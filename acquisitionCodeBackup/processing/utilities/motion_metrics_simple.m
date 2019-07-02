
cy = zeros(1,size(Y,4));
sz = size(Y);
Ys = reshape(Y,sz(1)*sz(2)*sz(3),sz(4));
M = mean(Ys,2);

for t=1:sz(4)
    disp(t)
    cy(t) = corr(M,Ys(:,t));
end