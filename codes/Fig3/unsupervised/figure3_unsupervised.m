 
clear;
clc;
N=16; % Node Number
totalnum=N*(N-1)/2;
outlier_ratio=0;
ii=[];
jj=[];
for k=2:N
    jj = [jj;ones(k-1,1)*k];
    ii = [ii;(1:(k-1))'];
end
edge = [ii,jj];
m = N*(N-1)/2;
d0 = sparse([1:m;1:m]',[ii;jj],[ones(1,m),-ones(1,m)],m,N);
for time=1:1000
    greedy_list = greedysampling(N);
    s=(0:N-1)'/(N-1);
    s_bar = s-mean(s);
    s_bar = s_bar(randperm(N));
    p0 = (s_bar(ii)-s_bar(jj)+1)/2;
    p0 = p0*(1-outlier_ratio) + (1-p0)*outlier_ratio;
    p1 = (s_bar(greedy_list(:,1))-s_bar(greedy_list(:,2))+1)/2;
    p1 = p1*(1-outlier_ratio) + (1-p1)*outlier_ratio;
    index=20:120;
    alist=[];
    zlist=zeros(totalnum,1);
    random_sample=[]; %random sampling 
    active_sample=[]; %unsupervised active sampling
    random_score=zeros(N,1);
    random_without_score=zeros(N,1);
    active_score=zeros(N,1);
    
    b=randperm(totalnum);
    for sample_i=1:totalnum
        p_re=rand(1);
        if p_re > p1(sample_i)
            active_sample=[active_sample
                greedy_list(sample_i,2:-1:1)];
        else
            active_sample=[active_sample
                greedy_list(sample_i,:)];
        end
        
        a=ceil(totalnum*rand(1));
        alist=[alist,a];
        p_re=rand(1);
        if p_re > p0(a)
            random_sample=[random_sample
                edge(a,2:-1:1)];
            zlist(sample_i)=1;
        else
            random_sample=[random_sample
                edge(a,:)];
        end
        
    end
    
    i = [1:totalnum,1:totalnum];
    j = [random_sample(:,1);random_sample(:,2)];
    k = [ones(1,totalnum),-ones(1,totalnum)];
    d1 = full(sparse(i,j,k,totalnum,N));
    
    j = [active_sample(:,1);active_sample(:,2)];
    d3= full(sparse(i,j,k,totalnum,N));
    
    for i = 1:length(index)
        sample_i = index(i);
        [random_score,flag] = lsqr(d1(1:sample_i,:)'*d1(1:sample_i,:),d1(1:sample_i,:)'*ones(sample_i,1));
        [active_score,flag] = lsqr(d3(1:sample_i,:)'*d3(1:sample_i,:),d3(1:sample_i,:)'*ones(sample_i,1));
        
        random_score = random_score-mean(random_score);
        active_score = active_score-mean(active_score);
        
        % kendall distance calculation
        Kendall_random(i)=corr(random_score,s_bar,'type','Kendall');
        Kendall_active(i)=corr(active_score,s_bar,'type','Kendall');
    end
    
    random_total(time,:)=Kendall_random;
    active_total(time,:)=Kendall_active;
end

for number=1:length(index)
    
    sd_greedy(number) = std(active_total(:,number));
    mean_greedy(number)=mean(active_total(:,number));
    
    sd_random(number) = std(random_total(:,number));
    mean_random(number)=mean(random_total(:,number));
    
end

ratio_temp= index;
plot(ratio_temp', mean_random','b','LineWidth',2 );hold on;
plot(ratio_temp',mean_greedy','r','LineWidth',2 );hold on;
legend('Random sampling','Unsupervised active sampling','fontsize',24);
xlabel('Sample number'); ylabel('Kendall');





