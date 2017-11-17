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

para=struct('reg', 0.01, 'maxiter', 100, 's0', 0, 'verbose', true, 'tol', 1e-4);

opt_s=struct('Method', 'lbfgs', 'DISPLAY', 0, 'MaxIter', 300, 'optTol', 1e-5, 'progTol', 1e-7);
base_s=minFunc(@func_s, s_init, opt_s,  ones(n_anno,1), para, pair);
base_kendall=calc_kendall(doc_diff, base_s, eps);
plot(base_s, doc_diff,  'b*');

s_init=randn(n_obj,1);
[s,alpha, obj, iter]=alter(s_init, alpha_init, pair, para);
kendall=calc_kendall(doc_diff, s, eps);
p=exp(s)/sum(exp(s));
plot(s, doc_diff,  'r.');
%kendall=corr(doc_diff, p, 'type', 'Kendall');