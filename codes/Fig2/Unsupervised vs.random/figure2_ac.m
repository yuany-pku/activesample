clc
clear

data=[];
N=16;
totalpair=N*(N-1)/2;
for i=1:N-1
    for j=i+1:N
        data=[data
            [i j]];
    end
end

time_number=100;
beta0_random=zeros(time_number, totalpair);
beta1_random=zeros(time_number, totalpair);
beta0_random2=zeros(time_number, totalpair);
beta1_random2=zeros(time_number, totalpair);
for time=1:time_number
    osting_list=greedysampling(N);
    random_sample=[];
    random_sample2=[];
    bb=randperm(totalpair);
    
    for sample_i=1:totalpair
        a=randperm(totalpair);
        random=data(a,:);
        random_sample=[random_sample
            random(sample_i,:)];
        random_sample2=[random_sample2
            osting_list(sample_i,:)];
        [beta0_random(time,sample_i) beta1_random(time,sample_i), beta2_random(time,sample_i)]=comp_betti(random_sample,N);
        [beta0_random2(time,sample_i) beta1_random2(time,sample_i),beta2_random2(time,sample_i)]=comp_betti(random_sample2,N);
        
    end
end

for number=1:totalpair
    mean_beta0_random(number)=mean(beta0_random(:,number));
    mean_beta0_random2(number)=mean(beta0_random2(:,number));
    mean_beta1_random(number)=mean(beta1_random(:,number));
    mean_beta1_random2(number)=mean(beta1_random2(:,number));
end

subplot(2,1,1)
plot(1:totalpair,mean_beta0_random2,'r-o',1:totalpair,mean_beta1_random2,'-b*','LineWidth',2);
legend('Betti0','Betti1');
title('Unsupervised active sampling');

subplot(2,1,2)
plot(1:totalpair,mean_beta0_random,'r-o',1:totalpair,mean_beta1_random,'-b*','LineWidth',2);
legend('Betti0','Betti1');
title('Random sampling');
 
