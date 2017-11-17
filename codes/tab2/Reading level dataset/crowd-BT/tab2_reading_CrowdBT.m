
clc;
clear;
tic
load('readdataset.mat')
load('groundtruthdata.mat')

n_obj = length(groundtruthdata);
run_time = 100;
totalnum = 4000;
s_bar=groundtruthdata-mean(groundtruthdata);
n_data=size(data_ref,1);

data_cx = [ones(n_data,1),data_ref(:,2:3)];

%% Test
active_para=struct('c', 0.1, 'kappa', 1e-4, 'taylor', 0, 'gamma',5, 'calc_iter', 1, ...
    'sel_method', 'greedy', 'anno_threshold', 1e-4,'max_iter', totalnum);

Tau = zeros(run_time,totalnum);
Acc = zeros(run_time,totalnum);
for time=1:run_time
    mu=zeros(n_obj,1);
    sigma=ones(n_obj,1);
    alpha=10;
    beta=1;
    tic()
    [mu, sigma, alpha, beta, auc, hist]...
        =active_learning(data_cx, mu, sigma, alpha, beta, s_bar, active_para);
    toc()
    Tau(time,:) = hist.tau';
    Acc(time,:) = hist.acc';
end
mean_Tau = mean(Tau)
mean_Acc = mean(Acc)
save('Result.mat');
toc

