function changeBasisDct=dct(patchSize)
% INPUT:
%   patchSize - size of the patches.
%
% OUTPUT:
%   changeBasisDct - patchSize^2 x patchSize^2-1 change-of-basis matrix. If
%       points is a matrix of patches in the pixel basis, then we multiply
%       by changeBasisDct to switch to the DCT basis:
%           pointsDct=points*changeBasisDct
%   
% NOTES:
%   See "The Nonlinear Statistics of High-Contrast Patches in Natural 
%   Images" by A.B. Lee, K.S. Pedersen, and D. Mumford [LPM03]. The matrices D, A,
%   and lambda below are the generalizations to arbitrary patch size of 
%   matrices D, A, and lambda in [LPM03]. I set changeBasisDct=A*lambda'. 
%   The change-of-basis equation given above is just the transpose of
%   equation 7 on page 89 of [LPM03].
%       v=lambda*A'*y
%       v'=y'*A*lambda'
%       v'=y'*changeBasisDct
%   The DCT basis vectors in matrix A are ordered differently in [LPM03]
%   then as produced by the code below.

D=Dgen(patchSize);
A=zeros(patchSize^2,patchSize^2-1);
v=zeros(patchSize^2,1);
for i=0:patchSize-1
    for j=0:patchSize-1
        for k=0:patchSize-1
            for m=0:patchSize-1
                v(1+k+patchSize*m,1)=cos(pi/patchSize*(k+1/2)*i)*cos(pi/patchSize*(m+1/2)*j);
            end
        end
        if i~=0 || j~=0
            A(:,i+patchSize*j)=v/sqrt(v'*D*v);
        end
    end
end
lambda=sparse(patchSize^2-1,patchSize^2-1);
for i=1:patchSize^2-1
   lambda(i,i)=1/norm(A(:,i))^2; 
end
changeBasisDct=A*lambda';



function D=Dgen(patchSize)
% INPUT:
%   patchSize - size of the patches.
%
% OUTPUT:
%   D - patchSize^2 x patchSize^2 contrast norm matrix. The 3 x 3 version
%   is given on page 88 of [LPM03]

D=sparse(patchSize^2,patchSize^2);
for i=1:patchSize^2 % loop over rows of D
    count=0;
    if i>patchSize
        D(i,i-patchSize)=-1;
        count=count+1;
    end
    if mod(i,patchSize)~=1
        D(i,i-1)=-1;
        count=count+1;
    end
    if mod(i,patchSize)~=0
        D(i,i+1)=-1;
        count=count+1;
    end
    if i<=patchSize^2-patchSize
        D(i,i+patchSize)=-1;
        count=count+1;
    end
    D(i,i)=count;    
end