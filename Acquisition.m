
function [IndexScenario] = Acquisition(Fn,problem,options)
f=problem.fun;
betaHat = problem.HBeta; 
switch options
%% --- select a scenario maximizes the similarity 
    case 'ASIM'% continuous design space; grid design space
        design  = problem.scenario; 
        [n,dim] = size(design);       
        beta1   = zeros(n,dim+1); beta0 = zeros(n,dim+1); 
    parfor i=1:n
        beta1(i,:)  = MLE([Fn;1,design(i,:),1],problem.HBeta,'EqWeight'); % m-by-1 vector
        beta0(i,:)  = MLE([Fn;1,design(i,:),0],problem.HBeta,'EqWeight'); 
    end 
        diff = max(abs([beta1(:,2:end);beta0(:,2:end)]-problem.HBeta(2:end)'),[],2);
        s = max(prctile(diff,5),0.05); 
        [x1,x2,x3]  = ndgrid(-1:s:1,-1:s:1,-1:s:1);
        pro = problem; pro.DS = [ x1(:) x2(:) x3(:)];
        %% calculate similarity
        Sim = zeros(n,1);
        X = [ones(n,1),design];
        p = problem.fun(X',problem.HBeta);   
        parfor i=1:n
           Sim(i) = p(i)*SimMeasure(problem.HBeta,beta1(i,:)',pro) + (1-p(i))*SimMeasure(problem.HBeta,beta0(i,:)',pro);   
        end
        [~,IndexScenario] = min(Sim);
        [N,~]=size(pro.DS);
%% --- select a scenario maximizes the similarity 
    case 'SMC'
        % N=100; X = -1 + 2*rand(N,1); 
        d=length(problem.TBeta)-1; q = 10^7; Q=10^5; Xq=[];sq=0;
        while sq<Q
            set = -1 + 2*rand(q,d);
            Pr = problem.fun([ones(q,1),set]',problem.HBeta);
            e =problem.eta+0.1>=Pr & Pr>=problem.eta-0.1;  Xq= [Xq;set(e,:)];  
            [sq,~]=size(Xq);
            if sq==0
                Xq=set(1:Q,:); sq=Q;
            end
        end
        Xq= Xq(1:Q,:); pro =problem;        pro.DS=Xq;
        
        lb = -ones(1,d); ub = -lb; 
        obj=@(x) eva(x,Fn,pro);
        options = optimoptions('surrogateopt','Display','off','PlotFcn',[],'MaxFunctionEvaluations',600);
        IndexScenario = surrogateopt(obj,lb,ub,options);
       
    case 'SML'
        design =problem.scenario; 
        Sim = zeros(size(design,1),1);
        Pr = problem.fun([ones(size(problem.DS,1),1),problem.DS]',problem.HBeta);
        e =problem.eta+0.1>=Pr & Pr>=problem.eta-0.1; pro =problem; 
        pro.DS=problem.DS(e,:);
    parfor i=1:size(design,1)
        Sample1       = [Fn;1,design(i,:),1];
        beta1  = MLE(Sample1,problem.HBeta,'EqWeight'); % m-by-1 vector
        Sample0       = [Fn;1,design(i,:),0];
        beta0  = MLE(Sample0,problem.HBeta,'EqWeight'); 
        p          = problem.fun([1,design(i,:)]',problem.HBeta);
        Sim(i)    = p*SimMeasure(problem.HBeta,beta1,pro) + (1-p)*SimMeasure(problem.HBeta,beta0,pro);  
    end
        [~,IndexScenario] = min(Sim);
        [N,~]=size(pro.DS);
    case 'SMs' % small-scale design space
        design =problem.scenario; pro =problem; 
        Sim = zeros(size(design,1),1);
    parfor i=1:size(design,1)
        Sample1       = [Fn;1,design(i,:),1];
        beta1  = MLE(Sample1,problem.HBeta,'EqWeight'); % m-by-1 vector
        Sample0       = [Fn;1,design(i,:),0];
        beta0  = MLE(Sample0,problem.HBeta,'EqWeight'); 
        p          = problem.fun([1,design(i,:)]',problem.HBeta);
        Sim(i)    = p*SimMeasure(problem.HBeta,beta1,pro) + (1-p)*SimMeasure(problem.HBeta,beta0,pro);  
    end
        [~,IndexScenario] = min(Sim);
        
%% ---A-optimality sampling: choose the point with smallest contribution to A-criterion
    case 'A-opt' 
    scenario = problem.scenario;
    KGquantity = zeros(size(scenario,1),1);
    parfor i=1:size(scenario,1)
        X = [1,scenario(i,:)];
        p = f(X',problem.HBeta);
        W = diag(p.*(1-p));
        KGquantity(i,1) = trace(X'*W*X);
    end
    [~,IndexScenario]=max(KGquantity);

%% --- D-optimality minimizes  the determinant of the information matrix
    case 'D-opt' 
    scenario = problem.scenario;
    KGquantity = zeros(size(scenario,1),1);
    parfor i=1:size(scenario,1)
        X = [1,scenario(i,:)];
        p = f(X',problem.HBeta);
        W = diag(p.*(1-p));
        KGquantity(i,1) = det(X'*W*X);
    end
    [~,IndexScenario]=max(KGquantity);

%% --- Randomly sample scenario
    case 'Random'       
    IndexScenario = randi(size(problem.scenario,1));
    
    case 'Uncertain'
    KGquantity = zeros(size(problem.scenario,1),1);
    for i=1:size(problem.scenario,1)
        KGquantity(i,1) = abs(f([1,problem.scenario(i,:)]',betaHat)-0.5);
    end
    [~,IndexScenario]=min(KGquantity); 
    
    case 'UTS' % sample the points nearest to the threshold value
    KGquantity = zeros(size(problem.scenario,1),1);
    for i=1:size(problem.scenario,1)
        KGquantity(i,1) = abs(f([1,problem.scenario(i,:)]',betaHat)-problem.eta);
    end
    [~,IndexScenario]=min(KGquantity);   
%% Expected Gradient Length: query the instance that would impart the greatest ...    
    case 'EGL'
    KGquantity = zeros(size(problem.scenario,1),1);
    parfor i=1:size(problem.scenario,1)
        x = [1,problem.scenario(i,:)]';
        KGquantity(i,1)=f(x,betaHat)*norm((1-f(x,betaHat))*x) +...
            (1-f(x,betaHat))*norm((0-f(x,betaHat))*x);
    end
    [~,IndexScenario]=max(KGquantity);   % maximum impact 
    case 'SMh' 
        Q = zeros(size(problem.scenario,1),1);
        for i=1:size(problem.scenario,1)
            Q(i,1) = abs(f([1,problem.scenario(i,:)]',betaHat)-problem.eta);
        end
        Qt = Q(Q<=0.2); ind=randi(length(Qt));
        IndexScenario = find(Qt(ind)==Q);
        if length(IndexScenario)>1
            IndexScenario=IndexScenario(randi(length(IndexScenario)));
        end
end
end


