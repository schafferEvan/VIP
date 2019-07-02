function registerMatfileDualColor(fullfilename)
addpath(genpath( '/rigel/theory/users/ess2129/code/NoRMCorre/'))
savepath = '/rigel/theory/users/ess2129/registered/';
parpool('local',12)

% parse 'fullfilename' into filename and path
u = strfind(fullfilename,'/');
filepath = fullfilename(1:u(end));
filename = fullfilename(u(end)+1:end);
m = matfile([filepath,filename]);
landmarkThreshold = 1000;

% load the data *** (colors flipped)
R = m.Y;
Y = m.R;
sizY = size(R);
mask = rand(sizY(1),sizY(2),sizY(3)) < 0.5; 
%mask = (max(R,[],4)<2).*(rand(sizY(1),sizY(2),sizY(3))<0.5); 
templates = m.templates;
ntemplates = size(m,'templates'); ntemplates = ntemplates(4);
save([savepath,filename(1:end-4),'reg.mat'],'R','Y','templates','landmarkThreshold','-v7.3'); %,'shiftStruct'
mfile = matfile([savepath,filename(1:end-4),'reg.mat'],'Writable',true);
clear m
movefile([filepath,filename],[filepath,'done/',filename])
disp(filename)



for j=1:ntemplates
    
    % select template
    template = templates(:,:,:,j).*(templates(:,:,:,j)>landmarkThreshold) + mask;    
    
    % coarse correction
    options = NoRMCorreSetParms('d1',sizY(1),'d2',sizY(2),'d3',sizY(3)','grid_size',[48,48,75],'bin_width',10,'mot_uf',2,'us_fac',10,...
    'method',{'median','mean'},'overlap_pre',24,'overlap_post',12,'max_dev',[10,10,10],...
    'mem_batch_size',10,'use_parallel',false,...
    'init_batch',10,'iter',1,'max_shift',50,'shifts_method','cubic','upd_template',false);

    %if ntemplates<3; options.iter=2; end

    Rth = R.*(R>landmarkThreshold);
    Rth = Rth+repmat(mask,1,1,1,sizY(4));
    [Rth,shifts] = normcorre_batch_even(Rth,options,template);
    
    Y = apply_shifts(Y,shifts,options);
    mfile.Y = Y;
    mfile.options = options;
    %clear Y
    
    R = apply_shifts(R,shifts,options);
    mfile.R = R;
    %clear R
    
    
    
    % fine correction
    options = NoRMCorreSetParms('d1',sizY(1),'d2',sizY(2),'d3',sizY(3)','grid_size',[20,36,20],'bin_width',10,'mot_uf',2,'us_fac',10,...
        'method',{'median','mean'},'overlap_pre',12,'overlap_post',12,'max_dev',[5,5,5],...
        'mem_batch_size',10,'use_parallel',false,...
        'init_batch',10,'iter',1,'max_shift',10,'shifts_method','cubic','upd_template',false);
    
    template = templates(:,:,:,j) + mask;
    Rth = R;  %.*(R>landmarkThreshold);
    Rth = Rth+repmat(mask,1,1,1,sizY(4));
    [Rth,shifts] = normcorre_batch_even(Rth,options,template);
    
    Y = apply_shifts(Y,shifts,options);
    mfile.Y = Y;
    mfile.options = options;
    %clear Y
    
    R = apply_shifts(R,shifts,options);
    mfile.R = R;
    %clear R

end



