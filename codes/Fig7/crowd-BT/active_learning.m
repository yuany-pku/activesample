function [mu, sigma, alpha, beta, auc, hist]...
    =active_learning(data, mu, sigma, alpha, beta, true_score, para)

gamma=getOpt(para,'gamma', 0);
calc_iter=getOpt(para,'calc_iter', 100);
anno_threshold=getOpt(para,'anno_threshold', 1e-4);
verbose=getOpt(para,'verbose', true);
sel_method=getOpt(para,'sel_method', 'greedy');
max_iter=getOpt(para,'max_iter', 120);

n_obj=length(mu);
data_temp=data; %找出里面还有的distinct pairs
Data=data;
Z = zeros(n_obj,n_obj);
for i=1:length(data_temp)
    Z(data_temp(i,2),data_temp(i,3)) = Z(data_temp(i,2),data_temp(i,3))+1;
    a = min(data_temp(i,2:3));
    b = max(data_temp(i,2:3));
    data_temp(i,2)= a;
    data_temp(i,3)= b;
end
edge=unique(data_temp,'rows');
w = zeros(1,length(edge));
p0 = w;

for i=1:length(edge)
    w(i) = Z(edge(i,2),edge(i,3)) + Z(edge(i,3),edge(i,2));
    p0(i) = Z(edge(i,2),edge(i,3))/w(i);
end

data = edge;
n_data=size(data,1);
auc=zeros(n_data,1);

[score, try_result]=init_score(data, mu, sigma, alpha, beta, para);

hist=struct('seq', zeros(n_obj,1), 'score', ones(n_obj,1),'mu',zeros(n_data,n_obj),'tau',zeros(max_iter,1),'acc',zeros(max_iter,1));
J_true = (true_score*ones(1,n_obj) - ones(n_obj,1)*true_score')>0;
num_J = sum(sum(J_true));

candidate=true(n_data,1);
for iter=1:max_iter
    if mod(iter,100)==0
        iter
    end
    %% select the highest score --- greedy algorithm
    if strcmp(sel_method, 'greedy')
        [hist.score(iter)]=max(score);
        idx=find(score==hist.score(iter));
        if length(idx)==1
            r=idx;
        else
            r=idx(randsample(length(idx),1));
        end
    elseif strcmp(sel_method, 'multinomial')
        r=find(mnrnd(1,score./sum(score)));
    elseif strcmp(sel_method, 'random')
        r=randsample(find(candidate),1);
    end
    
    hist.seq(iter)=r;
    %candidate(r)=false;
    %score(r)=0;
    
    %%  reveal the results i>_k j or i<_k j  and update parameters
    i=data(r, 2);
    j=data(r, 3);
    k=data(r, 1);
    if (rand(1)<p0(r))
        mu(i)=try_result{r,1}.mu1;
        mu(j)=try_result{r,1}.mu2;
        sigma(i)=try_result{r,1}.sigma1;
        sigma(j)=try_result{r,1}.sigma2;
        if abs(try_result{r,1}.alpha-alpha(k))<anno_threshold && abs(try_result{r,1}.beta-beta(k))<anno_threshold
            update_list=find( ((data(:,2)==i) | (data(:,3)==j) | (data(:,2)==j) | (data(:,3)==i)) & candidate);
        else
            update_list=find( (data(:,1) ==k | (data(:,2)==i) | (data(:,3)==j)| (data(:,2)==j) | (data(:,3)==i)) & candidate);
        end
        alpha(k)=try_result{r,1}.alpha;
        beta(k)=try_result{r,1}.beta;
    else
        mu(i)=try_result{r,2}.mu2;
        mu(j)=try_result{r,2}.mu1;
        sigma(i)=try_result{r,2}.sigma2;
        sigma(j)=try_result{r,2}.sigma1;
        if abs(try_result{r,2}.alpha-alpha(k))<anno_threshold && abs(try_result{r,2}.beta-beta(k))<anno_threshold
            update_list=find( ((data(:,2)==i) | (data(:,3)==j) | (data(:,2)==j) | (data(:,3)==i)) & candidate);
        else
            update_list=find( (data(:,1) ==k | (data(:,2)==i) | (data(:,3)==j)| (data(:,2)==j) | (data(:,3)==i)) & candidate);
        end
        alpha(k)=try_result{r,2}.alpha;
        beta(k)=try_result{r,2}.beta;
    end
    hist.mu(iter,:) = mu;
    
    if mod(iter,calc_iter)==0
        hist.tau(iter) = corr(mu,true_score,'type','Kendall');
        J = (mu*ones(1,n_obj) - ones(n_obj,1)*mu')>0;
        %hist.acc(iter) = calc_auc(true_score, mu);
        hist.acc(iter) = sum(sum(J&J_true))/num_J;
    end
    
    %% Update the new score and try_result
    for rr= 1:length(update_list)
        r=update_list(rr);
        i=data(r,2);
        j=data(r,3);
        k=data(r,1);
        [try_result{r,1}.mu1, try_result{r,1}.mu2, try_result{r,1}.sigma1, try_result{r,1}.simga2, try_result{r,1}.alpha,  try_result{r,1}.beta,...
            KL_win_o, KL_win_a, win_prob]=online_update(mu(i), mu(j), sigma(i), sigma(j), alpha(k), beta(k), para);
        [try_result{r,2}.mu1, try_result{r,2}.mu2, try_result{r,2}.sigma1, try_result{r,2}.simga2, try_result{r,2}.alpha,  try_result{r,2}.beta,...
            KL_lose_o, KL_lose_a, lose_prob]=online_update(mu(j), mu(i), sigma(j), sigma(i), alpha(k), beta(k), para);
        score(r)=win_prob*(KL_win_o+gamma*KL_win_a)+lose_prob*(KL_lose_o+gamma*KL_lose_a);
    end
end

end