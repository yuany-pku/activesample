
clc;
clear;
tic
N=16; % Node Number

run_time = 100;
totalnum = 120;

for ref=1:15
    str=strcat('.\PC-IQA dataset\','data',num2str(ref),'.mat');
    load(str); % load reference
    ground_truth_score=Hodgerank(data_ref);
    s_bar=ground_truth_score-mean(ground_truth_score);
    n_data=size(data_ref,1);
    data_cx = [ones(n_data,1),data_ref];
    n_anno=max(data_cx(:,1));
    n_obj=max(max(data_cx(:,2:3)));
    
    
    %% Test
    active_para=struct('c', 0.1, 'kappa', 1e-4, 'taylor', 0, 'gamma',5, 'calc_iter', 1, ...
        'sel_method', 'greedy', 'anno_threshold', 1e-4,'max_iter', totalnum);
    
    Tau = zeros(run_time,totalnum);
    for time=1:run_time
        mu=zeros(n_obj,1);
        sigma=ones(n_obj,1);
        alpha0=10;
        beta0=1;
        alpha=alpha0.*ones(n_anno,1);
        beta=beta0.*ones(n_anno,1);
        
        [mu, sigma, alpha, beta, auc, hist]...
            =active_learning(data_cx, mu, sigma, alpha, beta, s_bar, active_para);
        Tau(time,:) = hist.tau';
    end
    mean_Tau = mean(Tau);
    str_result=strcat('ref',num2str(ref),'.mat');
    save(str_result);
end
toc