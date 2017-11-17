%% Read data
clear;
data=dlmread('./read_data/all_pair.txt');

anno_quality=dlmread('./read_data/annotator_info.txt');
anno_quality=anno_quality(:,3);
doc_diff=dlmread('./read_data/doc_info.txt');
doc_diff=doc_diff(:,2);

n_anno=max(data(:,1));
n_obj=max(max(data(:,2:3)));
n_data=size(data,1);
%% set up initial parametmers 

mu=zeros(n_obj,1);
sigma=ones(n_obj,1);
alpha0=10;
beta0=1;
alpha=alpha0.*ones(n_anno,1);
beta=beta0.*ones(n_anno,1);
idx=randperm(n_data);
auc=zeros(n_obj,1);

online_para=struct('c', 0.1, 'kappa', 1e-4, 'taylor', 0);
for r=1: n_data
    i=data(idx(r), 2);
    j=data(idx(r), 3);
    k=data(idx(r), 1);
    [mu(i), mu(j), sigma(i), sigma(j), alpha(k), beta(k)]...
        =online_update(mu(i), mu(j), sigma(i), sigma(j), alpha(k), beta(k), online_para);
    if mod(r,100)==0
        auc(r)=calc_auc(doc_diff, mu);
        fprintf('Iter=%d, AUC=%.4f\n', r, auc(r));
    end
end
save(['exp_result/online_result_a_', num2str(alpha0),'_b_', num2str(beta0),'.mat'], 'mu', 'sigma', 'alpha', 'beta', 'idx', 'auc');