function [cx,cy,cz,cg] = sampleObject(v, n, noPlanes)

addpath('../Reconstruction');
addpath('../Reconstruction/utils');

% Simulate acquisition of 2D images
planes = arbitrarySlices(noPlanes,v,n,'random');

% create matrices with coordinates and grey-level intensities from
% structure 'planes'
for idx = 1:noPlanes,
    cx(:,:,idx) = planes(idx).x;
    cy(:,:,idx) = planes(idx).y;
    cz(:,:,idx) = planes(idx).z;
    cg(:,:,idx) = planes(idx).g;
end

close all;

