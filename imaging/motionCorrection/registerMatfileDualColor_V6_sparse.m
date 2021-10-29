function registerMatfileDualColor_v6_sparse(fullfilename)
addpath(genpath( '/moto/axs/users/code/NoRMCorre/'))
savepath = '/moto/axs/projects/registered/';
pp=parpool('local',12)
pp.IdleTimeout = 120;

% parse 'fullfilename' into filename and path
u = strfind(fullfilename,'/');
filepath = fullfilename(1:u(end));
filename = fullfilename(u(end)+1:end);
m = matfile([filepath,filename]);
landmarkThreshold = 200;

% load the data *** (colors flipped)
R = m.Y;
Y = m.R;
sizY = size(R);
mask = rand(sizY(1),sizY(2),sizY(3)) < 0.5;
%mask = (max(R,[],4)<2).*(rand(sizY(1),sizY(2),sizY(3))<0.5);
templates = m.templates;

templates = templates(:,:,:,[1,end]);

ntemplates = size(templates,4); %size(m,'templates'); ntemplates = ntemplates(4);
save([savepath,filename(1:end-4),'reg.mat'],'R','Y','templates','landmarkThreshold','-v7.3'); %,'shiftStruct'
mfile = matfile([savepath,filename(1:end-4),'reg.mat'],'Writable',true);
clear m Y
try
movefile([filepath,filename],[filepath,'done/',filename])
catch
mkdir([filepath,'done/'])
movefile([filepath,filename],[filepath,'done/',filename])
end
disp(filename)

iters=4

for j=1:ntemplates

    % for linear step with sparse line, use crosstalk
    landmarkThreshold = 200; %500
    landmarkThreshold_low = 110; %500
    % select template
    template = templates(:,:,:,j).*(templates(:,:,:,j)>landmarkThreshold_low) + mask;
    template(template>landmarkThreshold)=landmarkThreshold;

    % linear correction
    d1=sizY(1); d2=sizY(2); d3=sizY(3);
    options = NoRMCorreSetParms('d1',d1,'d2',d2,'d3',d3,'grid_size',[d1,d2,d3],'bin_width',10,'mot_uf',2,'us_fac',10,...
    'method',{'median','mean'},'overlap_pre',[floor(d1/4),floor(d2/4),floor(d3/4)],'overlap_post',[floor(d1/4),floor(d2/4),floor(d3/4)],'max_dev',[d1/4,d2/4,d3/4],...
    'mem_batch_size',10,'use_parallel',false,...
    'init_batch',10,'iter',1,'max_shift',50,'shifts_method','cubic','upd_template',false);

    %if ntemplates<3; options.iter=2; end

    Rth = R.*(R>landmarkThreshold_low);
    Rth(Rth>landmarkThreshold)=landmarkThreshold;
    clear R
    Rth = Rth+repmat(mask,1,1,1,sizY(4));
    [Rth,shifts] = normcorre_batch_even(Rth,options,template);
    clear Rth

    Y = apply_shifts(mfile.Y,shifts,options);
    mfile.Y = Y;
    mfile.options = options;
    clear Y

    R = apply_shifts(mfile.R,shifts,options);
    mfile.R = R;
    %clear R
    clear shifts


    landmarkThreshold = 400; %500
    landmarkThreshold_low = 200; %500
    % select template
    template = templates(:,:,:,j).*(templates(:,:,:,j)>landmarkThreshold_low) + mask;
    template(template>landmarkThreshold)=landmarkThreshold;

    % % coarse correction
    options = NoRMCorreSetParms('d1',sizY(1),'d2',sizY(2),'d3',sizY(3)','grid_size',[floor(d1/3),floor(d2/3),floor(d3/3)],'bin_width',10,'mot_uf',2,'us_fac',10,...
        'method',{'median','mean'},'overlap_pre',[floor(d1/6),floor(d2/6),floor(d3/6)],'overlap_post',[floor(d1/6),floor(d2/6),floor(d3/6)],'max_dev',[5,5,5],...
        'mem_batch_size',10,'use_parallel',false,...
        'init_batch',10,'iter',2,'max_shift',50,'shifts_method','cubic','upd_template',false);

    %if ntemplates<3; options.iter=2; end

    Rth = R.*(R>landmarkThreshold_low);
    Rth(Rth>landmarkThreshold)=landmarkThreshold;
    clear R
    Rth = Rth+repmat(mask,1,1,1,sizY(4));
    [Rth,shifts] = normcorre_batch_even(Rth,options,template);
    clear Rth

    Y = apply_shifts(mfile.Y,shifts,options);
    mfile.Y = Y;
    mfile.options = options;
    clear Y

    R = apply_shifts(mfile.R,shifts,options);
    mfile.R = R;
    %clear R
end



landmarkThreshold = 400; %500
landmarkThreshold_low = 110; %500
% select template
template = templates(:,:,:,end).*(templates(:,:,:,end)>landmarkThreshold_low) - landmarkThreshold_low + mask;
template(template>landmarkThreshold)=landmarkThreshold;


for j=1:iters
    
    % %medium correction 1
    options = NoRMCorreSetParms('d1',sizY(1),'d2',sizY(2),'d3',sizY(3)','grid_size',[35,35,35],'bin_width',10,'mot_uf',2,'us_fac',10,...
        'method',{'median','mean'},'overlap_pre',28,'overlap_post',28,'max_dev',[5,5,5],...
        'mem_batch_size',10,'use_parallel',false,...
        'init_batch',10,'iter',1,'max_shift',10,'shifts_method','cubic','upd_template',false);

    Rth = R.*(R>landmarkThreshold_low) - landmarkThreshold_low;
    Rth(Rth>landmarkThreshold)=landmarkThreshold;
    clear R
    Rth = Rth+repmat(mask,1,1,1,sizY(4));
    [Rth,shifts] = normcorre_batch_even(Rth,options,template);
    clear Rth

    Y = apply_shifts(mfile.Y,shifts,options);
    mfile.Y = Y;
    mfile.options = options;
    clear Y

    R = apply_shifts(mfile.R,shifts,options);
    mfile.R = R;



    % %medium correction 2
    options = NoRMCorreSetParms('d1',sizY(1),'d2',sizY(2),'d3',sizY(3)','grid_size',[33,33,33],'bin_width',10,'mot_uf',2,'us_fac',10,...
        'method',{'median','mean'},'overlap_pre',25,'overlap_post',25,'max_dev',[3,3,3],...
        'mem_batch_size',10,'use_parallel',false,...
        'init_batch',10,'iter',1,'max_shift',10,'shifts_method','cubic','upd_template',false);

    Rth = R.*(R>landmarkThreshold_low) - landmarkThreshold_low;
    Rth(Rth>landmarkThreshold)=landmarkThreshold;
    clear R
    Rth = Rth+repmat(mask,1,1,1,sizY(4));
    [Rth,shifts] = normcorre_batch_even(Rth,options,template);
    clear Rth

    Y = apply_shifts(mfile.Y,shifts,options);
    mfile.Y = Y;
    mfile.options = options;
    clear Y

    R = apply_shifts(mfile.R,shifts,options);
    mfile.R = R;


    % %medium correction 3
    options = NoRMCorreSetParms('d1',sizY(1),'d2',sizY(2),'d3',sizY(3)','grid_size',[30,30,30],'bin_width',10,'mot_uf',2,'us_fac',10,...
        'method',{'median','mean'},'overlap_pre',23,'overlap_post',23,'max_dev',[3,3,3],...
        'mem_batch_size',10,'use_parallel',false,...
        'init_batch',10,'iter',1,'max_shift',10,'shifts_method','cubic','upd_template',false);

    Rth = R.*(R>landmarkThreshold_low) - landmarkThreshold_low;
    Rth(Rth>landmarkThreshold)=landmarkThreshold;
    clear R
    Rth = Rth+repmat(mask,1,1,1,sizY(4));
    [Rth,shifts] = normcorre_batch_even(Rth,options,template);
    clear Rth

    Y = apply_shifts(mfile.Y,shifts,options);
    mfile.Y = Y;
    mfile.options = options;
    clear Y

    R = apply_shifts(mfile.R,shifts,options);
    mfile.R = R;
end



















