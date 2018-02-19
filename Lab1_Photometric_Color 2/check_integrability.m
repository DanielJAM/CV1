function [ p, q, SE ] = check_integrability( normals )
%CHECK_INTEGRABILITY check the surface gradient is acceptable
%   normals: normal image
%   p : df / dx
%   q : df / dy
%   SE : Squared Errors of the 2 second derivatives

% initalization
p = zeros(size(normals));
q = zeros(size(normals));
SE = zeros(size(normals));

% ========================================================================
% YOUR CODE GOES HERE
% Compute p and q, where
% p measures value of df / dx
% q measures value of df / dy
size_normals = size(normals);
for ht = 1: size_normals(1);
    for wd = 1: size_normals(2);
        p(ht,wd,1) = normals(ht,wd,1)/normals(ht,wd,3);
        q(ht,wd,1) = normals(ht,wd,2)/normals(ht,wd,3);
    end
end


% ========================================================================



p(isnan(p)) = 0;
q(isnan(q)) = 0;



% ========================================================================
% YOUR CODE GOES HERE
% approximate second derivate by neighbor difference
% and compute the Squared Errors SE of the 2 second derivatives SE


% ========================================================================
dpdy = zeros(size_normals(1)-1,size_normals(2)-1,1) 
dqdx = zeros(size_normals(1)-1,size_normals(2)-1,1)
for ht = 1: size_normals(1)-1;
    for wd = 1: size_normals(2)-1;
        dpdy(ht,wd,1) = (normals(ht,wd+1,1)/normals(ht,wd+1,3))-(normals(ht,wd,1)/normals(ht,wd,3));
        dqdx(ht,wd,1) = (normals(ht+1,wd,2)/normals(ht+1,wd,3))-(normals(ht,wd,2)/normals(ht,wd,3))
    end
end
 
SE = (dpdy-dqdx).^2 



end

