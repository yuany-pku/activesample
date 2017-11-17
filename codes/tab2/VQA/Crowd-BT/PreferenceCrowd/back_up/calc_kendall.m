function kendall=calc_kendall(x, y, epsilon)

    if ~exist('epsilon','var')
        epsilon=0;
    end
    
    n=length(x);
    cnt=0;
    for i=1:n-1
        for j=i+1: n
            if (x(i)-x(j))*(y(i)-y(j))<-abs(epsilon)
                cnt=cnt+1;
            end
        end
    end
    
    kendall=cnt*2/(n*(n-1));
                
end 