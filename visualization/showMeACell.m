
function showMeACell(A,dims,idx)

Ar = reshape( full(max(A,[],2)), dims(1),dims(2),dims(3));
Ar1 = max(Ar,[],3);
mx = max(Ar1(:));
R = Ar1*.7/mx;
G = Ar1*.7/mx;
B = Ar1*.7/mx;

Ar = reshape( full(max(A(:,idx),[],2)), dims(1),dims(2),dims(3));
Ar1 = max(Ar,[],3);
G = G - Ar1*.7/mx; % make idx cell red
B = B - Ar1*.7/mx;

mx = max([max(R(:)),max(G(:)),max(B(:))]);
R = R/mx; G = G/mx; B = B/mx;
RGB = cat(3,R,G,B);
figure; imagesc(RGB)