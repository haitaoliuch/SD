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

end







