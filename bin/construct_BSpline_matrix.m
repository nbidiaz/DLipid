function A = construct_BSpline_matrix(T,p,t0,te,dr)
% 
% Dynamics of the ethanolamine glycerophospholipid lipid remodeling network
%
% Copyright Mar 2012, Lu Zhang, all rights reserved.
% Boston College

    n = length(T);
    h = (te-t0)/p;
    S = (-1:p+2)*h+t0;
    A = zeros(n,p+4);
    for i = 1:n
        for j = -1:p+2
            A(i,j+2) = B(T(i),S(j+2),h,dr);
        end
    end
    
    A = [A(:,2)+2*A(:,1) A(:,3)-A(:,1) A(:,4:p) A(:,p+1)-A(:,p+3) A(:,p+2)+2*A(:,p+3)];

return
