%除以np实现归一化 
load result_eig_betti_new.mat
for time=1:20
    for i=1:120
        lamda2_active_total(time,i)=lamda2_active_total(time,i)/(16*i/120);
        lamda2_random_total(time,i)=lamda2_random_total(time,i)/(16*i/120);
    end
end

for number=1:120
    
        quan_osting=quantile(lamda2_active_total(:,number),[0.25,0.5,0.75]);
        quan_random=quantile(lamda2_random_total(:,number),[0.25,0.5,0.75]);

        median_osting(number)=quan_osting(2);
        error_L_osting(number)=quan_osting(2)-quan_osting(1);
        error_H_osting(number)=quan_osting(3)-quan_osting(2);

        median_random(number)=quan_random(2);
        error_L_random(number)=quan_random(2)-quan_random(1);
        error_H_random(number)=quan_random(3)-quan_random(2);

     
end


ratio_temp=1:number;
errorbar(ratio_temp',median_osting',error_L_osting',error_H_osting','r','LineWidth',1 );legend('Active sampling','Random sampling','fontsize',16);  hold on;
errorbar(ratio_temp',median_random',error_L_random',error_H_random','b','LineWidth',1 );legend('Active sampling','Random sampling','fontsize',16); hold on;
 %%%%betti errorbar