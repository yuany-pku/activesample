clear;
addpath(genpath('./minFunc_2012'));
data=dlmread('./read_data/all_pair.txt');
anno_quality=dlmread('./read_data/annotator_info.txt');
anno_quality=anno_quality(:,3);
doc_diff=dlmread('./read_data/doc_info.txt');
doc_diff=doc_diff(:,2);

n_anno=max(data(:,1));
n_obj=max(max(data(:,2:3)));
pair=cell(n_anno,1);
for i=1:n_anno
    pair{i}=data(data(:,1)==i, 2:3);
end

s_init=zeros(n_obj,1);
alpha_init=ones(n_anno,1);

para=struct('reg_0', 1, 'reg_s', 0, 'reg_alpha', 0,  'maxiter', 100, 's0', 0,...
             'uni_weight', true, 'verbose', true, 'tol', 1e-5);

opt_s=struct('Method', 'lbfgs', 'DISPLAY', 0, 'MaxIter', 300, 'optTol', 1e-5, 'progTol', 1e-7);
base_s=minFunc(@func_s, s_init, opt_s,  ones(n_anno,1), para, pair);
base_auc=calc_auc(doc_diff, base_s)
%base_kendall=calc_kendall(doc_diff, base_s, eps);
%plot(base_s, doc_diff,  'b*');

s_init=rand(n_obj,1);
[s,alpha, obj, iter]=alter(s_init, alpha_init, pair, para);
auc=calc_auc(doc_diff, s)

%kendall=calc_kendall(doc_diff, s, eps);
%plot(1:length(p), doc_diff(sort_idx), '*')
%kendall=corr(doc_diff, p, 'type', 'Kendall');
%plot(s, doc_diff,  'r.');