function [a b c]= comp_betti(incomp,N)

 
video_num=N;  
edge=incomp; 
triangle=[]; 
eps=1e-4;
sigma=1;

[temp1,temp2]=size(incomp);

Z = zeros(N,N); % paired comparison matrix
for k = 1:temp1
    a = incomp(k,:);
    Z(a(1),a(2)) = Z(a(1),a(2)) + 1;
end

thresh=0;
G=((Z+Z')>thresh);

for i=1:video_num,
    for j=(i+1):video_num,
        for k=(j+1):video_num,
            if ((G(i,j)>0)&&(G(j,k)>0)&&(G(k,i)>0)),
                triangle = [triangle
                    [i,j,k]];
            end
        end
    end
end


startJPlex

numEdges=size(edge,1);
numTriangles=size(triangle,1);

s2=ExplicitStream;
s2.add([1:video_num]',[0:video_num-1]');
s2.add(edge(:,:),[video_num:video_num-1+numEdges]');
if numTriangles>0
    s2.add(triangle(:,:),[video_num+numEdges:video_num-1+numEdges+numTriangles]');
end
s2.close

int2=Plex.Persistence.computeIntervals(s2);
betti_number=Plex.FilterInfinite(int2);
temp3=char(toString(betti_number));
s=sprintf('%c',temp3);
idx=strfind(s,'{');
idx2=strfind(s,'}');
temp3=s(idx+1:idx2-1);
betti_number=str2num(temp3);
if length(betti_number)==1
    betti_number(2)=0;
end
a=betti_number(1); % number of connected components
b=betti_number(2); % number of loops
c=numTriangles; % number of triangles
 