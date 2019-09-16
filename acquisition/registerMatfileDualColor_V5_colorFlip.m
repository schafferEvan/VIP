function registerMatfileDualColor_v3(fullfilename)
addpath(genpath( '/moto/axs/users/code/NoRMCorre/'))
savepath = '/moto/axs/projects/registered/';
pp=parpool('local',12)
pp.IdleTimeout = 120;

% parse 'fullfilename' into filename and path
u = strfind(fullfilename,'/');
filepath = fullfilename(1:u(end));
filename = fullfilename(u(end)+1:end);
m = matfile([filepath,filename]);

% load the data *** (colors flipped)
R = m.Y;
Y = m.R;
templates = m.templates;
templates = templates(:,:,:,[1,end]);


%disp(' ********* cropping ********** ')
%R = R(:,28:254,:,:);
%Y = Y(:,28:254,:,:);
%templates = templates(:,28:254,:,:);

sizY = size(R);
mask = rand(sizY(1),sizY(2),sizY(3)) < 0.5;
ntemplates = size(templates,4); %size(m,'templates'); ntemplates = ntemplates(4);
save([savepath,filename(1:end-4),'reg.mat'],'R','Y','templates','-v7.3'); %,'shiftStruct'
mfile = matfile([savepath,filename(1:end-4),'reg.mat'],'Writable',true);
clear m R
try
movefile([filepath,filename],[filepath,'done/',filename])
catch
mkdir([filepath,'done/'])
movefile([filepath,filename],[filepath,'done/',filename])
end
disp(filename)



for j=1:ntemplates

    landmarkThreshold = 400;
    % select template
    template = templates(:,:,:,j).*(templates(:,:,:,j)>landmarkThreshold) + mask;

    % linear correction
    d1=sizY(1); d2=sizY(2); d3=sizY(3);
    options = NoRMCorreSetParms('d1',d1,'d2',d2,'d3',d3,'grid_size',[d1,d2,d3],'bin_width',10,'mot_uf',2,'us_fac',10,...
    'method',{'median','mean'},'overlap_pre',[floor(d1/4),floor(d2/4),floor(d3/4)],'overlap_post',[floor(d1/4),floor(d2/4),floor(d3/4)],'max_dev',[d1/4,d2/4,d3/4],...
    'mem_batch_size',10,'use_parallel',false,...
    'init_batch',10,'iter',1,'max_shift',50,'shifts_method','cubic','upd_template',false);

    %if ntemplates<3; options.iter=2; end

    Rth = Y.*(Y>landmarkThreshold);
    clear Y
    Rth = Rth+repmat(mask,1,1,1,sizY(4));
    [Rth,shifts] = normcorre_batch_even(Rth,options,template);
    clear Rth

    R = apply_shifts(mfile.R,shifts,options);
    mfile.R = R;
    clear R
    
    Y = apply_shifts(mfile.Y,shifts,options);
    mfile.Y = Y;
    mfile.options = options;    
    clear shifts

    landmarkThreshold = 300;
    % select template
    template = templates(:,:,:,j).*(templates(:,:,:,j)>landmarkThreshold) + mask;

    % coarse correction
    options = NoRMCorreSetParms('d1',sizY(1),'d2',sizY(2),'d3',sizY(3)','grid_size',[floor(d1/3),floor(d2/3),floor(d3/3)],'bin_width',10,'mot_uf',2,'us_fac',10,...
    'method',{'median','mean'},'overlap_pre',[floor(d1/6),floor(d2/6),floor(d3/6)],'overlap_post',[floor(d1/6),floor(d2/6),floor(d3/6)],'max_dev',[15,15,15],...
    'mem_batch_size',10,'use_parallel',false,...
    'init_batch',10,'iter',2,'max_shift',50,'shifts_method','cubic','upd_template',false);

    %if ntemplates<3; options.iter=2; end

    Rth = Y.*(Y>landmarkThreshold);
    clear Y
    Rth = Rth+repmat(mask,1,1,1,sizY(4));
    [Rth,shifts] = normcorre_batch_even(Rth,options,template);
    clear Rth

    R = apply_shifts(mfile.R,shifts,options);
    mfile.R = R;
    clear R
    
    Y = apply_shifts(mfile.Y,shifts,options);
    mfile.Y = Y;
    mfile.options = options;    
    clear shifts


    landmarkThreshold = 250;
    % select template
    template = templates(:,:,:,j).*(templates(:,:,:,j)>landmarkThreshold) + mask;


    % medium correction
    options = NoRMCorreSetParms('d1',sizY(1),'d2',sizY(2),'d3',sizY(3)','grid_size',[25,36,25],'bin_width',10,'mot_uf',2,'us_fac',10,...
        'method',{'median','mean'},'overlap_pre',12,'overlap_post',12,'max_dev',[10,10,10],...
        'mem_batch_size',10,'use_parallel',false,...
        'init_batch',10,'iter',2,'max_shift',10,'shifts_method','cubic','upd_template',false);

    
    Rth = Y.*(Y>landmarkThreshold);
    clear Y
    Rth = Rth+repmat(mask,1,1,1,sizY(4));
    [Rth,shifts] = normcorre_batch_even(Rth,options,template);
    clear Rth

    R = apply_shifts(mfile.R,shifts,options);
    mfile.R = R;
    clear R
    
    Y = apply_shifts(mfile.Y,shifts,options);
    mfile.Y = Y;
    mfile.options = options;    
    clear shifts

    landmarkThreshold = 200;
    % select template
    template = templates(:,:,:,j).*(templates(:,:,:,j)>landmarkThreshold) + mask;


    % finest correction
    options = NoRMCorreSetParms('d1',sizY(1),'d2',sizY(2),'d3',sizY(3)','grid_size',[20,30,20],'bin_width',10,'mot_uf',2,'us_fac',10,...
        'method',{'median','mean'},'overlap_pre',8,'overlap_post',8,'max_dev',[5,5,5],...
        'mem_batch_size',10,'use_parallel',false,...
        'init_batch',10,'iter',1,'max_shift',10,'shifts_method','cubic','upd_template',false);

    
    Rth = Y.*(Y>landmarkThreshold);
    clear Y
    Rth = Rth+repmat(mask,1,1,1,sizY(4));
    [Rth,shifts] = normcorre_batch_even(Rth,options,template);
    clear Rth

    R = apply_shifts(mfile.R,shifts,options);
    mfile.R = R;
    clear R
    
    Y = apply_shifts(mfile.Y,shifts,options);
    mfile.Y = Y;
    mfile.options = options;    
    clear shifts

end

