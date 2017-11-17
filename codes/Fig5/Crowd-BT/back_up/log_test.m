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
cnt=0;
for i=1:n_anno
    pair{i}=data(data(:,1)==i, 2:3);
    cnt=cnt+length(pair{i});
end

pos_J=[];
neg_J=[];
for i=1:n_anno
    pos_J=[pos_J; pair{i}(:,1)];
    neg_J=[neg_J; pair{i}(:,2)];    
end

X=sparse([1:cnt,1:cnt], [pos_J; neg_J], [ones(cnt,1); -ones(cnt,1)]);
y=-ones(cnt,1);
w_init = zeros(n_obj,1);
funObj = @(w)LogisticLoss(w,X,y);
options.Method='lbfgs';
[s,funVal]=minFunc(@LogisticLoss,w_init,options,X,y);
