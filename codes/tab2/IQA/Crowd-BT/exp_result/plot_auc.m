load active_greedy_g_1_a_5.mat
greedy_1_auc=auc;

load active_greedy_g_10.mat
greedy_10_auc=auc;

load active_random_g_10.mat
random_auc=auc;

load active_multinomial_g_10.mat
multinomial_auc=auc;

load online_result_a_10_b_1.mat
online_auc=auc;

figure
hold on
plot(find(greedy_1_auc), greedy_1_auc(greedy_1_auc>0), 'r-');
%plot(find(multinomial_auc), multinomial_auc(multinomial_auc>0), 'm-');
%plot(find(random_auc), random_auc(random_auc>0), 'b-');
plot(find(online_auc), online_auc(online_auc>0), 'k--');
hold off