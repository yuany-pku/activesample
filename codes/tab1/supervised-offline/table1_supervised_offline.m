clear;
clc;
warning off;
n=[16,20,24,28,32,100];
count=0;
for index=1:1:6
    N=n(index);
    count=count+1;
    initialnum = 0;
    totalnum=N*(N-1)/2;
    run_time =100;
    lambda = 1e-1;
    outlier_ratio=0;
    data_model = 1;
    solve_model = 1;
    % 1:uniform;  2:Bradley-Terry  3:Thurstone-Mosteller   4:Angular transform
    eps = 1e-4;
    active_Kendall_distance=zeros(run_time,totalnum);
    edge=[];
    for j=2:N
        for i=1:(j-1)
            edge=[edge
                i j];
        end
    end
    i = edge(:,1);
    j = edge(:,2);
    m = N*(N-1)/2;
    d = sparse([1:m;1:m]',[i;j],[ones(1,m),-ones(1,m)],m,N);
    a=0;
    for time=1:run_time
        edge_index = [];
        s=(0:N-1)'/(N-1);
        s_bar = s-mean(s);
        s_bar = s_bar(randperm(N));
        s = s_bar;
        p0 = score2prob(s, data_model);
        tic();
        w = zeros(1,m);
        pi = w;
        if (initialnum>0)
            temp=randperm(N*(N-1)/2);
            w(temp(1:initialnum)) = 1;
            pi(temp(1:initialnum)) = double(rand(initialnum,1) < p0(temp(1:initialnum)));
        end
        w_random=w;
        pi_random = pi;
        dis_active=[];
        [s_initial,Sigma0] = solveHodge(d,w,pi,solve_model,lambda);
        for num=1:totalnum
            p = score2prob(s_initial, solve_model);
            edge_informationgain=zeros(m,3);
            for ii=1:(N*(N-1)/2)
                w_new=w;
                w_new(ii) = w_new(ii) + 1;
                pi_new1 = pi;
                pi_new2 = pi;
                pi_new1(ii) = (pi(ii)*w(ii)+1)/(w(ii)+1);
                pi_new2(ii) = pi(ii)*w(ii)/(w(ii)+1);
                [s_new1,Sigma1] = solveHodge(d,w_new,pi_new1,solve_model,lambda);
                [s_new2,TTT] = solveHodge(d,w_new,pi_new2,solve_model,lambda);
                gain=p(ii)*(s_new1 - s_initial)'*Sigma0*(s_new1 - s_initial)+...
                    (1-p(ii))*(s_new2 - s_initial)'*Sigma0*(s_new2 - s_initial)+...
                    trace(Sigma0*pinv(full(Sigma1))) + log(det(Sigma1));
                edge_informationgain(ii,:)=[edge(ii,:) gain];
            end
            informationgain=edge_informationgain(:,3);
            [infor index]=sort(informationgain,'descend');
            maxgain_index=index(1);
            edge_index = [edge_index,index(1)];
            temp_p=rand(1);
            p_re=rand(1);
            if (temp_p<p0(maxgain_index) && p_re>outlier_ratio) || (temp_p> p0(maxgain_index) && p_re<outlier_ratio)
                pi(maxgain_index) = (pi(maxgain_index)*w(maxgain_index)+1)/(w(maxgain_index)+1);
            else
                pi(maxgain_index) = pi(maxgain_index)*w(maxgain_index)/(w(maxgain_index)+1);
            end
            w(maxgain_index) = w(maxgain_index)+1;
            [active_score,Sigma0] = solveHodge(d,w,pi,solve_model,lambda);
            s_initial = active_score;
            end
        a=a+toc();
    end
    A(count)=a;
end