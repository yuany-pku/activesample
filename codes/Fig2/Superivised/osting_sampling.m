clc;
clear;

load data6.mat %load ref1
load osting_list.mat %load osting sampling list

temp=randperm(16);
for i=1:3840
    for j=1:2
        data_ref(i,j)= temp(data_ref(i,j));
    end
end

%%%%%%%%%%%%%ground-truth score
ground_truth_score=Hodgerank(data_ref);
%%%%%%%%%%%%%%%%%%%%%%%%%
kendall_random=[];
kendall_osting=[];

for time=1:20
    a=randperm(3840);

    random=[];
    incomp=[];
    count=0;
    index=randperm(32);
    % kendall_osting=zeros(3840,1);
    for kk=1:32
        k=index(kk);
        for i=1:120
            for j=1:120
                if ((data_ref((k-1)*120+j,1)==osting_list(i,1) && data_ref((k-1)*120+j,2)==osting_list(i,2)) || (data_ref((k-1)*120+j,1)==osting_list(i,2) && data_ref((k-1)*120+j,2)==osting_list(i,1)))
                    incomp=[incomp
                        data_ref((k-1)*120+j,:)];

                end
            end
             if count>=150
            break;
             else
            osting_score=Hodgerank(incomp);

            count1=size(incomp);
            count=count1(1);
           
                kendall_osting(time,count)=corrkendall(ground_truth_score,osting_score);
                random=[random
                    data_ref(a(count),:)];
                random_score=Hodgerank(random);
                kendall_random(time,count)=corrkendall(ground_truth_score,random_score);
             end

        end
       

    end
   
end


for number=1:count
    for time=1:20
        kendall_osting_1(time)=kendall_osting(time,number);
        kendall_random_1(time)=kendall_random(time,number);
        quan_osting=quantile(kendall_osting_1,[0.25,0.5,0.75]);
        quan_random=quantile(kendall_random_1,[0.25,0.5,0.75]);

        median_osting(number)=quan_osting(2);
        error_L_osting(number)=quan_osting(2)-quan_osting(1);
        error_H_osting(number)=quan_osting(3)-quan_osting(2);

        median_random(number)=quan_random(2);
        error_L_random(number)=quan_random(2)-quan_random(1);
        error_H_random(number)=quan_random(3)-quan_random(2);

    end
end


ratio_temp=1:number;
errorbar(ratio_temp',median_osting',error_L_osting',error_H_osting','r','LineWidth',1 );legend('Active sampling','Random sampling','fontsize',16);  hold on;
errorbar(ratio_temp',median_random',error_L_random',error_H_random','b','LineWidth',1 );legend('Active sampling','Random sampling','fontsize',16);
save('result6.mat','median_osting','error_L_osting','error_H_osting','median_random','error_L_random','error_H_random','ratio_temp')