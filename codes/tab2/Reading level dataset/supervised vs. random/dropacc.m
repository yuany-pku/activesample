% DESCRIPTION:
%   This is the code to reproduce the result in the Fig.4.

tic;
clc;
clear;

N=491; % Node Number
initialnum = 0;

run_time =1;
lambda = 1e-1;
load readdataset.mat; % load reference
load groundtruthdata.mat;
s_bar=groundtruthdata(:,2);
J0 = (s_bar*ones(1,N)-ones(N,1)*s_bar') > 0;
totalnum=4000;
data_model =1;
solve_model = 1; %注意下面一并需要改
% 1:uniform;  2:Bradley-Terry  3:Thurstone-Mosteller   4:Angular transform
eps = 1e-4;
ref_num=1;
edge=[];
data_ref_temp=data_ref; %找出里面还有的distinct pairs
Z = zeros(N,N);
for i=1:length(data_ref)
    Z(data_ref(i,1),data_ref(i,2)) = Z(data_ref(i,1),data_ref(i,2))+1;
    a=min(data_ref(i,:));
    b=max(data_ref(i,:));
    data_ref_temp(i,1)=a;
    data_ref_temp(i,2)=b;
end
edge=unique(data_ref_temp,'rows');

i = edge(:,1);
j = edge(:,2);
ii = (i-1)*N + i;
jj = (j-1)*N + j;
ij = (j-1)*N + i;

i0=[];j0=[];
for k=2:N
    j0 = [j0;ones(k-1,1)*k];
    i0 = [i0;(1:(k-1))'];
end
edge_complete=[i0 j0];
[ccc,i_index]=ismember(edge, edge_complete,'rows');
p0=[];
for k = 1:length(i_index)
    p0 = [p0,Z(edge(k,1),edge(k,2))/(Z(edge(k,1),edge(k,2))+Z(edge(k,2),edge(k,1)))];
end

m = length(edge);
d = sparse([1:m;1:m]',[i;j],[ones(1,m),-ones(1,m)],m,N);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%  Input matrix with paired comparison, the matrix has
%  2 columns, and for each row, the rank of the first column
%  is higher than the second column for this comparison.
for time=1:run_time
    aaaaaaaaaaaaaaaaaaaaaaa=time
    edge_index = [];
    
    w = zeros(m,1);
    pi = w;
    
    w_random=w;
    pi_random = pi;
    edges_add=[];
    dis_random=[];
    dis_active=[];
    
    Sigma_0_inv=pinv(d'*((w*ones(1,N)).*d) + lambda*eye(N));
    s_initial = Sigma_0_inv * (d' * (w .* (2*pi-1)));
    s_initial = s_initial - mean(s_initial);
    s_random = s_initial;
    Sigma_0_inv_random = Sigma_0_inv;
    active_sample=[];
    active_match_ratio=zeros(totalnum,1);
    random_match_ratio=zeros(totalnum,1);
    active_match_ratio2=zeros(totalnum,1);
    random_match_ratio2=zeros(totalnum,1);
    for num=1:totalnum
        if (round(num/100) ==num/100 )
            display(num);
        end
        count0=0;
        count1=0;
        %             bbbbbbbbbbbb = num
        %%算每条边的pij.
        p = score2prob(s_initial, solve_model);
        p = p(i_index);
        
        d_S0_inv_d = Sigma_0_inv(ii) + Sigma_0_inv(jj) - 2*Sigma_0_inv(ij);  %sum(d*Sigma_0_inv.*d,2);
        d_s0 = d*s_initial;
        d_S1_inv_d = d_S0_inv_d - d_S0_inv_d.^2./(1+d_S0_inv_d);
        
        ds_Sigma0_ds = ((1 + d_s0).^2-4*d_s0.*p).*d_S0_inv_d./(1+d_S0_inv_d).^2;
        informationgain = ds_Sigma0_ds - log(1- d_S1_inv_d) - d_S1_inv_d;
        
        index=find(informationgain == max(informationgain));%找出information gain 最大的边
        maxgain_index = index(1);
        edge_index = [edge_index,maxgain_index];
        
        edges_add=[edges_add
            edge(maxgain_index,:)];
        
        
        %%%%%%%%%%%%%%%%% 加入某条边
        tep=rand(1);
        Sigma0_inv_d = Sigma_0_inv*d(maxgain_index,:)';
        if tep<p0(maxgain_index)
            active_score = s_initial + Sigma0_inv_d*(1-d_s0(maxgain_index))/(1+d_S0_inv_d(maxgain_index));
        else
            active_score = s_initial + Sigma0_inv_d*(-1-d_s0(maxgain_index))/(1+d_S0_inv_d(maxgain_index));
        end
        
        Sigma_0_inv = Sigma_0_inv - Sigma0_inv_d*Sigma0_inv_d'/(1+d_S0_inv_d(maxgain_index));
        s_initial = active_score;
        %算match ratio
        %         for jjj=1:len
        %             if active_score(data(jjj,1))- active_score(data(jjj,2))>0
        %                 count0=count0+1;
        %             end
        %         end
        %         active_match_ratio(num,1)=count0/len;
        %%%算kendall
        active_match_ratio2(num,1)=corr(active_score,s_bar,'type','Kendall');% l2 distance
        total_cnt=0;
        cnt=0;
        %           %%%%%%%%算ACC
        
        J1 = (active_score*ones(1,N)-ones(N,1)*active_score') > 0;
        active_match_ratio(num,1)=sum(sum(J0.*J1))/sum(sum(J0));
        %%%随机加入random
        random_temp=ceil(m*rand(1));
        tep = rand(1);
        
        Sigma_0_inv_random_d = Sigma_0_inv_random*d(random_temp,:)';
        if tep < p0(random_temp)
            random_score = s_random + Sigma_0_inv_random_d*(1-d(random_temp,:)*s_random)/(1+d(random_temp,:)*Sigma_0_inv_random_d);
        else
            random_score = s_random + Sigma_0_inv_random_d*(-1-d(random_temp,:)*s_random)/(1+d(random_temp,:)*Sigma_0_inv_random_d);
        end
        
        s_random = random_score;
        Sigma_0_inv_random = Sigma_0_inv_random - Sigma_0_inv_random_d * Sigma_0_inv_random_d'/(1+d(random_temp,:)*Sigma_0_inv_random_d);
        
        %         for jjj=1:len
        %             if random_score(data(jjj,1))- random_score(data(jjj,2))>0
        %                 count1=count1+1;
        %             end
        %         end
        %         random_match_ratio(num,1)=count1/len;
%         random_match_ratio2(num,1)=corr(random_score,s_bar,'type','Kendall');% l2 distance
%         
%         total_cnt=0;
%         cnt=0;
%         %%%%%%%%算ACC
%         J2 = (random_score*ones(1,N)-ones(N,1)*random_score') > 0;
%         random_match_ratio(num,1)=sum(sum(J0.*J2))/sum(sum(J0));;
    end
    
%     active2_match_ratio(time,:)=active_match_ratio';
%     random2_match_ratio(time,:)=random_match_ratio';
%     active2_match_ratio2(time,:)=active_match_ratio2';
%     random2_match_ratio2(time,:)=random_match_ratio2';
    %
end


% for number=1:totalnum
%     %     sd_random(number) = std(random_l2_distance(:,number));
%     mean_acc_active(number)=mean(active2_match_ratio(:,number));
%     mean_acc_random(number)=mean(random2_match_ratio(:,number));
%     mean_acc_active2(number)=mean(active2_match_ratio2(:,number));
%     mean_acc_random2(number)=mean(random2_match_ratio2(:,number));
%     
% end
% subplot(2,1,1)
% plot(20:totalnum', mean_acc_random(20:end)','b','LineWidth',2 );hold on;
% plot(20:totalnum',mean_acc_active(20:end)','r','LineWidth',2 );hold on;
% xlabel('Sample number');
% ylabel('ACC');
% legend('Random sampling','Active sampling','fontsize',24);
% subplot(2,1,2)
% plot(20:totalnum', mean_acc_random2(20:end)','b','LineWidth',2 );hold on;
% plot(20:totalnum',mean_acc_active2(20:end)','r','LineWidth',2 );hold on;
% xlabel('Sample number');
% ylabel('Kendall');
% legend('Random sampling','Active sampling','fontsize',24);
toc;
