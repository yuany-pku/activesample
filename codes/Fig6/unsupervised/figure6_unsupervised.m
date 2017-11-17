

clc;
clear;

data_ind=[];
N=16;
totalnum=N*(N-1)/2;
index_initial=floor((120*log(N))/N);
ref_num=2;
for i=1:N-1
    for j=i+1:N
        data_ind=[data_ind
            [j i]];
    end
end

for ref=1:ref_num;
    str=strcat('.\PC-IQA dataset\','data',num2str(ref),'.mat');
    load(str); % load reference 1
    %  Input matrix with paired comparison, the matrix has
    %  2 columns, and for each row, the rank of the first column
    %  is higher than the second column for this comparison.
    ground_truth_score=Hodgerank(data_ref);
    s_bar=ground_truth_score-mean(ground_truth_score);
    
    for time=1:100
        greedy_list = greedysampling(N);
        greedy_perm=randperm(N);
        for ii=1:totalnum
            greedy_list(ii,1)=greedy_perm(greedy_list(ii,1));
            greedy_list(ii,2)=greedy_perm(greedy_list(ii,2));
            
        end
        
        len=length(data_ref);
        a=randperm(totalnum);
        active_sample=[];
        random_sample=[];
        random_sample2=[];
        
        Kendall_active=zeros(1,totalnum);
        Kendall_random=zeros(1,totalnum);
        Kendall_random2=zeros(1,totalnum);
        
        for sample_i=1:totalnum
            active_score=zeros(N,1);
            random_score=zeros(N,1);
            random_without_score=zeros(N,1);
            
            b=randperm(totalnum);
            tep=randperm(len);
            data_ref=data_ref(tep,:);
            for i=1:len
                data=data_ref;
                if data(i,1)==data_ind(a(sample_i),1) && data(i,2)==data_ind(a(sample_i),2) || data(i,1)==data_ind(a(sample_i),2)&& data(i,2)==data_ind(a(sample_i),1)
                    random_sample=[random_sample
                        data(i,:)];
                    break;
                end
            end
            
            for i=1:len
                data=data_ref;
                if data(i,1)==data_ind(b(sample_i),1) && data(i,2)==data_ind(b(sample_i),2) || data(i,1)==data_ind(b(sample_i),2)&& data(i,2)==data_ind(b(sample_i),1)
                    random_sample2=[random_sample2
                        data(i,:)];
                    break;
                end
            end
            
            for i=1:len
                data=data_ref;
                if data(i,1)==greedy_list(sample_i,1) && data(i,2)==greedy_list(sample_i,2) || data(i,1)==greedy_list(sample_i,2) && data(i,2)==greedy_list(sample_i,1)
                    active_sample=[active_sample
                        data(i,:)];
                    break;
                end
            end
            
            
            if sample_i>=floor((totalnum*log(N))/N)
                active_score=Hodgerank(active_sample);
                active_score=active_score-mean(active_score);
            end
            
            if sample_i>=floor((totalnum*log(N))/N)
                random_without_score=Hodgerank(random_sample);
                random_without_score=random_without_score-mean(random_without_score);
            end
            if sample_i>=floor((totalnum*log(N))/N)
                random_score=Hodgerank(random_sample2);
                random_score=random_score-mean(random_score);
            end
            if length(random_without_score)==length(s_bar)&& length(s_bar)==length(active_score)&& length(random_score)==length(s_bar)
                Kendall_active(sample_i)=corr(active_score,s_bar,'type','Kendall');
                Kendall_random(sample_i)=corr(random_without_score,s_bar,'type','Kendall');
                Kendall_random2(sample_i)=corr(random_score,s_bar,'type','Kendall');
            end
            
            Kendall_random(1:index_initial)=NaN;
            Kendall_active(1:index_initial)=NaN;
            Kendall_random2(1:index_initial)=NaN;
        end
        
        active_total(time,:)= Kendall_active;
        random_total(time,:)= Kendall_random;
        random2_total(time,:)= Kendall_random2;
        
    end
    
    for number=1:totalnum
        mean_greedy(number)=mean(active_total(:,number));
        mean_random(number)=mean(random_total(:,number));
        mean_random2(number)=mean(random2_total(:,number));
    end
    
    ratio_temp= 1:number;
    ratio_temp=ratio_temp*N/(totalnum*log(N));
    str_result=strcat('ref',num2str(ref),'.mat');
    save(str_result);
end

