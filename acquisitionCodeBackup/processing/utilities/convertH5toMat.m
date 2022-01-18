
addpath(genpath('C:\Users\Axel-SCAPE\Documents\_code\analysis\normcorre\'))
expFolder = 'H:\2018_07_19_hungry\h5files\';%'G:\habaRegistered\spon\';


handle = 'reg*';
outerDir = dir([expFolder,handle]);
mkdir([expFolder,'mats'])


for k=1:length(outerDir)
    h5file = [expFolder,outerDir(k).name];
    Y = bigread2(h5file);
    if k==1
        template = Y(:,:,:,1);
    end
    save([expFolder,'mats\',outerDir(k).name(1:end-2),'mat'],'Y','template','-v7.3')
end
