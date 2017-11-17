% DESCRIPTION:
%    Greedy sampling for paired comparison data.

% INPUT ARGUMENTS:
% Num_Node    The number of node in the graph.
 
% OUTPUT ARGUMENTS:
% EDGES       The final sampling list. 
%
% DATE: 2014-5-15 
%
% REFERENCES: Enhanced statistical rankings via targeted data collection,ICML 2014
 
function EDGES = greedysampling(Num_Node)

n=Num_Node; N=n*(n-1)/2;

% define a minimally connected graph
w = zeros(N,1); w(1) = 1; kk = 1;
for ii = n-1:-1:2, kk = kk+ii; w(kk) = 1; end

% construct the gradient and a list of edges for the full laplacian
edges = 1:N; m = length(edges);
E = zeros(m,2);
ii = 1; % current row number
rr = 0; % number of elements in above rows
for kk = 1:m
    while edges(kk)>rr+(n-ii), rr = rr+(n-ii); ii = ii+1; end
    E(kk,1)=ii;
    E(kk,2)=edges(kk)-rr+ii; % column number
end
grad = sparse(1:m,E(:,2),ones(m,1),m,n)-sparse(1:m,E(:,1),ones(m,1),m,n);
lap = grad'*(sparse(1:N,1:N,w,N,N)*grad);

% now build the well-connected graph
m1 = nnz(w); %current number of edges
m2 = N; % desired number of edges
EDGES = zeros(m2,2); EDGES(1:m1,:) = E(w==1,:); % current edges
AC = zeros(m2-m1,1);

for ii = 1:m2-m1
    [V,D] = eigs(grad'*spdiags(w,0,N,N)*grad,[],2,'SA');
    vec2=V(:,2); AC(ii) = D(2,2);
    [XX,i] = sort(abs(grad*vec2),'descend');
    for jj = 1:N, % add the edge (i,j) which maximizes |v_i - v_j|
        if w(i(jj))==0, w(i(jj))=1; EDGES(m1+ii,:) = E(i(jj),:); break; end
    end
end
 
