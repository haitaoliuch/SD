%% for continuous design space
clc,clear,close all, warning off
for j=2:2
    if j==1
        method='SML' ;
    elseif j==2
        method='SMC' ;
    elseif j==3
        method='UTS' ;
    elseif j==4
        method='Uncertain' ;
    elseif j==5
        method='EGL';   
    elseif j==6
        method='D-opt' ;  
    elseif j==7
        method='A-opt' ;
    elseif j==8
        method='Random';
    end
    load('..\SLSim2\HPro\P_30D.mat')
    N =1000; R =1;    P= '30D';   pro.eta= 0.8; sim = zeros(R,N); %pro    = Problem(P);
    tic
for r=1:R
    problem=pro; [a,d]=size(problem.scenario); 
    problem.HBeta = zeros(d+1,1); Sample=[]; temp=[]; 
    for n=0:N-1 
        if n<d*2
           x_n = -1 + 2*rand(1,d);
        else  
           x_n = Acquisition(Sample,problem,method); 
        end  
       if unifrnd(0,1) <= problem.fun([1,x_n]',problem.TBeta) 
           Sample =[Sample;1,x_n,1];
       else
           Sample =[Sample;1,x_n,0];   
       end
        problem.HBeta = MLE(Sample,problem.HBeta,'eqweight');               % Update regression coefficient  based on observations         
        if mod(n+1,10)==0      
            temp(1,n+1) = SimMeasure(problem.TBeta,problem.HBeta, problem);  %---Measure Similarity  
        end
    end 
    sim(r,:)= temp;
end
toc
 Rep.SimMean  = mean(sim);  
folder= fullfile('..\SLSim2\Result',method);   [~,~] = mkdir(folder); % make new folder
save(fullfile(folder,sprintf('P_%s_N%d_E%d.mat',P,N,pro.eta*10)),'Rep','sim'); 
clearvars  -except j
end
