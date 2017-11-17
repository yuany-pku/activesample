% Compute the kl-divergence from the parameters of two dirichlet
% distributions.
%
% p - 1xN vector parameters for distibution 1
% q - 1xN vector parameters for distibution 2
% k - the resulting kl-divergence
function k=KLDirichlet(p, q)

k= sum((p-q).*(psi(p)-psi(sum(p)))) + sum(dGamma(q,p)) + dGamma(sum(p),sum(q));


function res=dGamma(p, q)
lis=[p;q];
[~, i]=min(lis,[],1);
res=zeros(1,size(lis,2));
for ik=1:size(lis,2)
    a=ceil(lis(i(ik),ik)/10)-1;
    pp=p(ik)-(a*10);
    qp=q(ik)-(a*10);
    inter=1;
    for k=1:a*10
        inter=inter*((p(ik)-k)/(q(ik)-k));
    end
    res(ik)= log(inter) + (gammaln(pp)- gammaln(qp));
end

