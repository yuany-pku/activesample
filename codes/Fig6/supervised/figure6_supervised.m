clc;
clear;

N=16; % Node Number
initialnum = 0;
totalnum=120 ;
run_time =100;
lambda = 1;
outlier_ratio=0;
solve_model =1;
% 1:uniform;  2:Bradley-Terry  3:Thurstone-Mosteller   4:Angular transform
eps = 1e-4;
ref_num=15;


for ref=1:ref_num; %reference index, there are 15 references in PC-IQA dataset
    str=strcat('.\PC-IQA dataset\','data',num2str(ref),'.mat');
    load(str); % load reference
    ground_truth_score=Hodgerank(data_ref);
    s_bar=ground_truth_score-mean(ground_truth_score);
    
    data_ref_temp=data_ref;
    Z = zeros(N,N);
    for i=1:length(data_ref)
        Z(data_ref(i,1),data_ref(i,2)) = Z(data_ref(i,1),data_ref(i,2))+1;
        data_ref_temp(i,1)=min(data_ref(i,:));
        data_ref_temp(i,2)=max(data_ref(i,:));
    end
    edge=unique(data_ref_temp,'rows');
    w = zeros(1,length(edge));
    p0 = w;
    for i=1:length(edge)
        w(i) = Z(edge(i,1),edge(i,2)) + Z(edge(i,2),edge(i,1));
        p0(i) = Z(edge(i,1),edge(i,2))/w(i);
    end
    
    i = edge(:,1);
    j = edge(:,2);
    ii = (i-1)*N + i;
    jj = (j-1)*N + j;
    ij = (j-1)*N + i;
    
    m = length(edge);
    d = sparse([1:m;1:m]',[i;j],[ones(1,m),-ones(1,m)],m,N);
    
    for time=1:run_time
        edge_index = [];
        Sigma_0_inv=1/lambda*eye(N);
        s_initial = zeros(N,1);
        s_initial = s_initial - mean(s_initial);
        s_random = s_initial;
        Sigma_0_inv_random = Sigma_0_inv;
        
        active_sample=[];
        active_match_ratio=zeros(totalnum,1);
        random_match_ratio=zeros(totalnum,1);
        for num=1:totalnum
            p = (s_initial(i) - s_initial(j) + 1)/2;
            d_S0_inv_d = Sigma_0_inv(ii) + Sigma_0_inv(jj) - 2*Sigma_0_inv(ij);
            d_s0 = d*s_initial;
            d_S1_inv_d = d_S0_inv_d - d_S0_inv_d.^2./(1+d_S0_inv_d);
            ds_Sigma0_ds = ((1 + d_s0).^2-4*d_s0.*p).*d_S0_inv_d./(1+d_S0_inv_d).^2;
            informationgain = ds_Sigma0_ds - log(1- d_S1_inv_d) - d_S1_inv_d;
            index=find(informationgain == max(informationgain));
            maxgain_index = index(1);
            edge_index = [edge_index,maxgain_index];
            
            tep=rand(1);
            Sigma0_inv_d = Sigma_0_inv*d(maxgain_index,:)';
            if tep<p0(maxgain_index)
                active_score = s_initial + Sigma0_inv_d*(1-d_s0(maxgain_index))/(1+d_S0_inv_d(maxgain_index));
            else
                active_score = s_initial + Sigma0_inv_d*(-1-d_s0(maxgain_index))/(1+d_S0_inv_d(maxgain_index));
            end
            
            Sigma_0_inv = Sigma_0_inv - Sigma0_inv_d*Sigma0_inv_d'/(1+d_S0_inv_d(maxgain_index));
            s_initial = active_score;
            
            active_Kendall(num,1)=corr(active_score,s_bar,'type','Kendall');
            random_temp = find(mnrnd(1,w/sum(w)));
            tep = rand(1);
            Sigma_0_inv_random_d = Sigma_0_inv_random*d(random_temp,:)';
            if tep < p0(random_temp)
                random_score = s_random + Sigma_0_inv_random_d*(1-d(random_temp,:)*s_random)/(1+d(random_temp,:)*Sigma_0_inv_random_d);
            else
                random_score = s_random + Sigma_0_inv_random_d*(-1-d(random_temp,:)*s_random)/(1+d(random_temp,:)*Sigma_0_inv_random_d);
            end
            
            s_random = random_score;
            Sigma_0_inv_random = Sigma_0_inv_random - Sigma_0_inv_random_d * Sigma_0_inv_random_d'/(1+d(random_temp,:)*Sigma_0_inv_random_d);
            random_Kendall(num,1)=corr(random_score,s_bar,'type','Kendall'); 
            
        end
        active_total(time,:)=active_Kendall';
        random_total(time,:)=random_Kendall';
    end
   
    mean_Kendall_active=mean(active_total);
    mean_Kendall_random=mean(random_total);

    str_result=strcat('ref',num2str(ref),'.mat');
    save(str_result);
end
