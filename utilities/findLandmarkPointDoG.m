
function [points,im_DoG] = findLandmarkPointDoG(im, params)

sz = size(im);
if nargin<2
    s1 = 0.8; %1;
    s2 = 1.8; %1.5;
    th = 40;
else
    if ~isfield(params,'th')
        params.th = 250;
    end
    th = params.th;
    s1 = params.s1;
    s2 = params.s2;
end


m1 = imgaussfilt3(im,s1);
m2 = imgaussfilt3(im,s2);

im_DoG = m1-m2;

validPts = find( (m1-m2 > th) );

[x,y,z] = ind2sub(sz, validPts );
points = [x,y,z];
