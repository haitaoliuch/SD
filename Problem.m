function pro = Problem(options)

%---Define the logistic regression function
pro.fun=@(x,beta) 1./(1+exp(-x' * beta));
switch options
    case '3D1S'
        pro.TBeta =[1;1;1;1];   
        step=0.1;    [x_3, x_2,x_1]=ndgrid(-1:step:1,-1:step:1,-1:step:1);
        pro.scenario =[ x_1(:) x_2(:) x_3(:)]; % collection of scenarios
        step = 0.01; [x_1, x_2,x_3] =ndgrid(-1:step:1,-1:step:1,-1:step:1);
        pro.DS =[ x_1(:) x_2(:) x_3(:)]; % collection of Designs
    case '3D2S'
        pro.TBeta =[1;1;1;1];  
        step=0.2; [x_3, x_2,x_1]=ndgrid(-1:step:1,-1:step:1,-1:step:1);
        pro.scenario =[ x_1(:) x_2(:) x_3(:)]; % collection of scenarios 
        step = 0.01; [x_1, x_2,x_3] =ndgrid(-1:step:1,-1:step:1,-1:step:1);
        pro.DS =[ x_1(:) x_2(:) x_3(:)]; % collection of Designs
    case 'S0'
        pro.TBeta =[1;1;1;1];  
        step=0.2; [x_3, x_2,x_1]=ndgrid(-1:step:1,-1:step:1,-1:step:1);
        pro.scenario =[ x_1(:) x_2(:) x_3(:)]; % collection of scenarios 
        step = 0.2; [x_1, x_2,x_3] =ndgrid(-1:step:1,-1:step:1,-1:step:1);
        pro.DS =[ x_1(:) x_2(:) x_3(:)]; % collection of Designs
    case 'S1'
        pro.TBeta =[1;1;1;1];  
        step=0.2; [x_3, x_2,x_1]=ndgrid(-1:step:1,-1:step:1,-1:step:1);
        pro.scenario =[ x_1(:) x_2(:) x_3(:)]; % collection of scenarios 
        step = 0.2; [x_1, x_2,x_3] =ndgrid(-1:step:1,-1:step:1,-1:step:1);
        pro.DS =[ x_1(:) x_2(:) x_3(:)]; % collection of Designs
        n=size(pro.DS,1);
        y = pro.fun([ones(n,1),pro.DS]',pro.TBeta);
        pro.eta = max(y)*0.8;
    case '4D2S'
        pro.TBeta =[1;1;1;1;1];  
        step=0.2; 
        [x_4,x_3, x_2,x_1]=ndgrid(-1:step:1,-1:step:1,-1:step:1,-1:step:1);
        pro.scenario =[ x_1(:) x_2(:) x_3(:) x_4(:)]; 
        
    case 'bank'      
        data = xlsread('bank.xlsx',1);     % import data
        pro.scenario =normalize((data(1:10000,1:end-1)),'range', [-1 1]);% Normalized Predictors
        pro.DS = pro.scenario;  pro.Y= data(1:10000,end);
        pro.TBeta = MLE([ones(size(pro.scenario,1),1) pro.DS pro.Y],[],'eqweight');
    case 'VCP'      
        data = xlsread('VCP.xlsx',2);     % import data
        pro.scenario =normalize((data(:,1:end-1)),'range', [-1 1]);% Normalized Predictors
        data = xlsread('VCP.xlsx',1);
        pro.DS = normalize((data(:,1:end-1)),'range', [-1 1]); [c,r]=size(pro.DS);
        pro.TBeta = ones(r+1,1);% MLE([ones(size(pro.DS,1),1) data],[],'eqweight');
        % pro.Y=data(:,end);  
        
    case 'MA2D'
        pro.TBeta =[1;1;1;0;0;0];   
        step=0.05;    [x_1, x_2]=ndgrid(-1:step:1,-1:step:1);
        pro.scenario =[ x_1(:) x_2(:) x_1(:).^2, x_2(:).^2, x_1(:).*x_2(:)]; % collection of scenarios
        step = 0.01; [x_1, x_2] =ndgrid(-1:step:1,-1:step:1);
        pro.DS =[ x_1(:) x_2(:) x_1(:).^2, x_2(:).^2, x_1(:).*x_2(:)]; % collection of Designs
        mu =0; sigma =2;
        pro.fun=@(x,beta) cdf('Normal',x' * beta,mu,sigma);  
    case '1D'
        pro.TBeta =[1;1];  
        step=0.1; [x_1]=ndgrid(-1:step:1);n=length(x_1); 
        for i=1:n
            pro.mu(i)=pro.fun([1,x_1(i)]',pro.TBeta);
        end 
        pro.scenario = x_1(:); % collection of scenarios 
        pro.DS =x_1(:); % collection of Designs  
    case '50D'
        d = 50; n = 10000; r=1;
        pro.TBeta =ones(d+1,1);  
        pro.scenario =randsphere(n,d,r); 
end







