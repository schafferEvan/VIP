

function template = makeManyMatfilesFromTiffStack(trialfolder,savePath,experimentName,chunkSize,zeroThreshold,timeToExclude,cropping,template)
if ismac
    addpath(genpath('/Users/evan/Dropbox/_code/NoRMCorre'))
else
end

if nargin<3
    chunkSize = 375;                                                        % number of time steps per matfile (375 steps * 8 parts = 3000)
    zeroThreshold = 110;                                                    % threshold for re-zeroing data
    timeToExclude = 300;                                                    % number of time steps at beginning to skip
    savePath = trialfolder;
end

d = dir([trialfolder,'*Fly*']);
u = strfind(d(1).name,'_');
handle = d(1).name(1:u(end));

nChunks = ceil((length(d)-timeToExclude)/chunkSize);
tStart = timeToExclude+1;
tEnd = tStart + chunkSize - 1;

overWrite = false;

for k=1:nChunks
    tRange = tStart:tEnd;
    filename = [experimentName,'_',handle,num2str(tRange(1)),'.mat'];
    
    if ~overWrite && exist([savePath,filename],'file')
        continue
    end
    
    Z = double(loadtiff([trialfolder,handle,num2str(tRange(1)),'.tiff']));
    if ~isempty(cropping.x); Z = Z(cropping.x,:,:); end
    if ~isempty(cropping.y); Z = Z(:,cropping.y,:); end
    if ~isempty(cropping.z); Z = Z(:,:,cropping.z); end
    sz = size(Z);
    Y = zeros(sz(1),sz(2),sz(3),length(tRange));
    Y(:,:,:,1)=Z;
    if (k==1)&&(nargin<8)
        template = Z;
    end
    for j=2:length(tRange)
        Z = double(loadtiff([trialfolder,handle,num2str(tRange(j)),'.tiff']));
        if ~isempty(cropping.x); Z = Z(cropping.x,:,:); end
        if ~isempty(cropping.y); Z = Z(:,cropping.y,:); end
        if ~isempty(cropping.z); Z = Z(:,:,cropping.z); end
        Y(:,:,:,j) = Z;
    end
    Y(Y<zeroThreshold)=zeroThreshold;
    Y(isnan(Y))=zeroThreshold; % update this in no0 function
    Y = Y-zeroThreshold;
    
    save([savePath,filename],'Y','template','cropping','-v7.3')
    
    tStart = tEnd+1;
    tEnd = min(tStart + chunkSize - 1,length(d));
end