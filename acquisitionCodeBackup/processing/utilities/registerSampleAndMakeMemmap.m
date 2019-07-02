
addpath(genpath('C:\Users\Axel-SCAPE\Documents\_code\analysis\normcorre\'))
expFolder = 'D:\2018_07_20_starved_54hr\';
outerFolder = [expFolder,'tiff_stacks\'];


handle = 'fly*';
outerDir = dir([outerFolder,handle]);


for k=1; %1:length(outerDir)
    datafolder = [outerFolder,outerDir(k).name,'\'];
    
    d = dir([datafolder,'*.tiff']);
    tRange = zeros(size(d));
    for j=1:length(d)
        u = strfind(d(j).name,'_');
        tRange(j) = str2double(d(j).name(u(end)+1:end-5));
        if isnan(tRange(j))
            tRange(j) = str2double(d(j).name(u(end-1)+1:u(end)-1));
        end
    end
    tRange = sort(tRange,'ascend');
    
    
    
    
    tRange = 401:500;
    
    
    
%     Z = double(loadtiff([datafolder,handle,num2str(tRange(1)),'.tiff']));
    Z = double(loadtiff([datafolder,d(1).name]));
    sz = size(Z);
    Y = zeros(sz(1),sz(2),sz(3),length(tRange));
    Y(:,:,:,1)=Z;
    for j=2:length(tRange)
        Y(:,:,:,j) = double(loadtiff([datafolder,d(j).name]));
    end
    Y(Y<110)=110;
    Y(isnan(Y))=110; % update this in no0 function
    Y = Y-110;
    
    sizY = size(Y);
    
    filename = outerDir(k).name;
    
    options = NoRMCorreSetParms('d1',sizY(1),'d2',sizY(2),'d3',sizY(3)','grid_size',[24,24,24],'bin_width',1,'mot_uf',2,'us_fac',10,...
        'method',{'median','mean'},'overlap_pre',24,'overlap_post',24,'max_dev',[20,20,20],...
        'output_type','h5','mem_batch_size',10,'h5_filename',[expFolder,'registered_',filename,'.h5'],'use_parallel',false,...
        'init_batch',1,'iter',1,'max_shift',50,'shifts_method','linear');
    
    save(filename,'Y','-v7.3')
    clear Y
    m = matfile(filename);
    
    %normcorre(Y,options);
    batch(@normcorre,0,{m,options})
end
