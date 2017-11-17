clear;
clc;
startJplex

% generate a graph with 16 convexes and all edges

edges=[];
for i = 1:16;
    for j = (i+1):16;
        edges=[edges;[i,j]];
    end
end

numEdges=size(edges,1);

% random sampling edges to add 

percent=0.7;
numSample=percent*numEdges;
a=randperm(numEdges);
samp=a(1:numSample);
incomp=edges(samp,:);

triangles=[];
for i=1:16,
    for j=(i+1):16,
        for k=(j+1):16,
            c=0;
            for p1=1:numSample,
                if incomp(p1,:)==[i,j],c=c+1;end
            end
            for p2=1:numSample,
                if incomp(p2,:)==[j,k],c=c+1;end
            end
            for p3=1:numSample,
                if incomp(p3,:)==[i,k],c=c+1;end
            end
            if c==3,
                triangles = [triangles;[i,j,k]];
            end
        end
    end
end

numTriangles=size(triangles,1);
% add edges and make graph for persistent homology
s2=ExplicitStream;
s2.add([1:16]',[0:15]');
s2.add(incomp(:,:),[16:15+numSample]');
%s2.ensure_all_faces
s2.add(triangles(:,:,:),[16+numSample:15+numSample+numTriangles]');
s2.close

int2=Plex.Persistence.computeIntervals(s2);
Plex.plot(int2,'Barcode',15+numSample+numTriangles);