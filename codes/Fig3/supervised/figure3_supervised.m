clear;
clc;
warning off;
tic;
N=16; % Node Number
initialnum = 0;
totalnum=N*(N-1)/2;
run_time = 1000;
lambda = 1e-1;
outlier_ratio=0;
data_model = 1;
solve_model = 1;
% 1:uniform;  2:Bradley-Terry  3:Thurstone-Mosteller   4:Angular transform
eps = 1e-4;

random_kendall_distance=zeros(run_time,totalnum);
active_kendall_distance=zeros(run_time,totalnum);
i=[];j=[];
for k=2:N
    j = [j;ones(k-1,1)*k];
    i = [i;(1:(k-1))'];
end
edge = [i,j];
m = N*(N-1)/2;
d = sparse([1:m;1:m]',[i;j],[ones(1,m),-ones(1,m)],m,N);

for time=1:run_time
    edge_index = [];
    s=(0:N-1)'/(N-1);
    s_bar = s-mean(s);
    s_bar = s_bar(randperm(N));
    s = s_bar;
    p0 = score2prob(s, data_model);
   
    
    w = zeros(m,1); 
    pi = w;  
    
    if (initialnum>0)
        temp=randperm(N*(N-1)/2);
        w(temp(1:initialnum)) = 1;
        pi(temp(1:initialnum)) = double(rand(initialnum,1) < p0(temp(1:initialnum)));
    end
    pi_random = pi;
    dis_random=[];
    dis_active=[];
    
    Sigma_0_inv=pinv(d'*((w*ones(1,N)).*d) + lambda*eye(N));
    s_initial = Sigma_0_inv * (d' * (w .* (2*pi-1)));
    s_initial = s_initial - mean(s_initial);
    s_random = s_initial;
    Sigma_0_inv_random = Sigma_0_inv;
    
    for num=1:totalnum
        p = score2prob(s_initial, solve_model); 
        ii = (i-1)*N + i;
        jj = (j-1)*N + j;
        ij = (j-1)*N + i;
        d_S0_inv_d = Sigma_0_inv(ii) + Sigma_0_inv(jj) - 2*Sigma_0_inv(ij);  
        d_s0 = d*s_initial;
        d_S1_inv_d = d_S0_inv_d - d_S0_inv_d.^2./(1+d_S0_inv_d);
        
        ds_Sigma0_ds = ((1 + d_s0).^2-4*d_s0.*p).*d_S0_inv_d./(1+d_S0_inv_d).^2;
        informationgain = ds_Sigma0_ds - log(1- d_S1_inv_d) - d_S1_inv_d;
        
        [infor index]=sort(informationgain,'descend');
        maxgain_index=index(1);  
        edge_index = [edge_index,index(1)];
         
        temp_p=rand(1);
        p_re=rand(1);
        Sigma0_inv_d = Sigma_0_inv*d(maxgain_index,:)';
        if (temp_p<p0(maxgain_index) && p_re>outlier_ratio) || (temp_p> p0(maxgain_index) && p_re<outlier_ratio)
            active_score = s_initial + Sigma0_inv_d*(1-d_s0(maxgain_index))/(1+d_S0_inv_d(maxgain_index));
        else
            active_score = s_initial + Sigma0_inv_d*(-1-d_s0(maxgain_index))/(1+d_S0_inv_d(maxgain_index));
        end
        Sigma_0_inv = Sigma_0_inv - Sigma0_inv_d*Sigma0_inv_d'/(1+d_S0_inv_d(maxgain_index));
        s_initial = active_score;
         
        random_temp=ceil(rand(1)*N*(N-1)/2);
        temp_p=rand(1);
        p_re=rand(1);
        Sigma_0_inv_random_d = Sigma_0_inv_random*d(random_temp,:)';
        if temp_p<p0(random_temp) && p_re>outlier_ratio || temp_p> p0(random_temp) && p_re<outlier_ratio
            random_score = s_random + Sigma_0_inv_random_d*(1-d(random_temp,:)*s_random)/(1+d(random_temp,:)*Sigma_0_inv_random_d);
        else
            random_score = s_random + Sigma_0_inv_random_d*(-1-d(random_temp,:)*s_random)/(1+d(random_temp,:)*Sigma_0_inv_random_d);
        end
        s_random = random_score;
        Sigma_0_inv_random = Sigma_0_inv_random - Sigma_0_inv_random_d * Sigma_0_inv_random_d'/(1+d(random_temp,:)*Sigma_0_inv_random_d);
        
   
        kendall_random = corr(random_score,s_bar, 'type','Kendall');
        dis_random=[dis_random
            kendall_random];
        kendall_active = corr(active_score,s_bar,'type', 'Kendall');
        dis_active=[dis_active
            kendall_active];
    end
    
    random_kendall_distance(time,:)=dis_random';
    active_kendall_distance(time,:)=dis_active';
 
end

for number=1:totalnum
    mean_random(number)=mean(random_kendall_distance(:,number));
    mean_active(number)=mean(active_kendall_distance(:,number));
end


plot(20:totalnum', mean_random(20:end)','b','LineWidth',2 );hold on;
plot(20:totalnum',mean_active(20:end)','r','LineWidth',2 );hold on;
xlabel('Sample number');
ylabel('kendall');
legend('Random sampling','Supervised active sampling','fontsize',24);
toc;
 