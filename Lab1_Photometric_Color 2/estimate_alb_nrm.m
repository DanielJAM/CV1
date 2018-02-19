function [ albedo, normal ] = estimate_alb_nrm( image_stack, scriptV, shadow_trick)
%COMPUTE_SURFACE_GRADIENT compute the gradient of the surface
%   image_stack : the images of the desired surface stacked up on the 3rd
%   dimension
%   scriptV : matrix V (in the algorithm) of source and camera information
%   shadow_trick: (true/false) whether or not to use shadow trick in solving
%   	linear equations
%   albedo : the surface albedo
%   normal : the surface normal


[h, w, ~] = size(image_stack);
if nargin == 2
    shadow_trick = true;
end

% create arrays for 
%   albedo (1 channel)
%   normal (3 channels)
albedo = zeros(h, w, 1);
normal = zeros(h, w, 3);

% =========================================================================
% YOUR CODE GOES HERE
% for each point in the image array
%   stack image values into a vector i
%   construct the diagonal matrix scriptI
%   solve scriptI * scriptV * g = scriptI * i to obtain g for this point
%   albedo at this point is |g|
%   normal at this point is g / |g|
    function [pointwise_albedo, pointwise_normal]= pointwise(height,width)
        size_stack = size(image_stack);
        i = image_stack(height, width, 1:size_stack(3));
        i=reshape(i, [1,size_stack(3)]);
        scriptI = eye(size_stack(3)).*i;
        if shadow_trick == true
            rhs = scriptI*i';
            lhs = scriptI* scriptV;
        else
            rhs = i';
            lhs = scriptV;
        end
        g=linsolve(lhs,rhs);
        pointwise_albedo = norm(g);
        pointwise_normal = g'./norm(g);
    end

for ht = 1:h
    for wd = 1:w
        albedo(ht,wd,1) = pointwise(ht,wd);
        normal(ht,wd,1:3) = pointwise(ht,wd);
    end
end


% =========================================================================

end

