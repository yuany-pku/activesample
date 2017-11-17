
clear;
clc;
data=[];
N=16;
greedy_list = greedysampling(N);
totalpair=N*(N-1)/2;
for i=1:N-1   
    data=[data;[i*ones(N-i,1),(i+1:N)']];
end

ind = floor((1:20)/20*totalpair);
iter=100;

lamda_with=zeros(iter,length(ind)); 
for time=1:iter
    a=ceil(totalpair*rand(max(ind),1));   
    i = [1:max(ind),1:max(ind)];
    j = [data(a,1);data(a,2)];
    k = [ones(1,max(ind)),-ones(1,max(ind))];
    d1 = sparse(i,j,k,max(ind),N);  
    for k = 1:length(ind)
        L = d1(1:ind(k),:)'*d1(1:ind(k),:);
        [V,D] = eigs(L,[],2,'SA');
        lamda1(k)=D(2,2);
    end
    lamda_with(time,:)=lamda1;
end
for number=1:length(ind)
    mean_lambda_with(number)=mean(lamda_with(:,number));
end

i = [1:max(ind),1:max(ind)];   
k = [ones(1,max(ind)),-ones(1,max(ind))];
j = [greedy_list(1:max(ind),1);greedy_list(1:max(ind),2)];
d_active = sparse(i,j,k,max(ind),N);
for k = 1:length(ind)
    L = d_active(1:ind(k),:)'*d_active(1:ind(k),:);
    [V,D] = eigs(L,[],2,'SA');
    lamda_active(k)=D(2,2);
    L = d_active(1:(ind(k)-1),:)'*d_active(1:(ind(k)-1),:);
    [V,D2] = eigs(L,[],2,'SA');
    delta = (d_active(ind(k),:)*V(:,2))^2;
    Dn = eigs(L,[],1,'LM');
    lamda_upper(k)=D2(2,2) + delta/(1+(2-delta)^2/(Dn-D2(2,2)));
end

plot(ind,mean_lambda_with,'g',ind,lamda_active,'b','LineWidth',2);


title('')
xlim([N,totalpair]);
ylim([0,N])
xlabel('Sample Number')
ylabel('\lambda_2')

legend('Random sampling','Unsupervised active sampling',...
    'fontsize',24);