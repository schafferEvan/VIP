function reregisterMatfileDualColor_V6(fullfilename)
addpath(genpath( '/moto/axs/users/code/NoRMCorre/'))
savepath = '/moto/axs/projects/registered/';
pp=parpool('local',12)
pp.IdleTimeout = 120;

% parse 'fullfilename' into filename and path
u = strfind(fullfilename,'/');
filepath = fullfilename(1:u(end));
filename = fullfilename(u(end)+1:end);
m = matfile([filepath,filename]);

% load the data *** (colors correct)
R = m.R;
Y = m.Y;
sizY = size(Y);
mask = rand(sizY(1),sizY(2),sizY(3)) < 0.5;
templates = m.templates;
templates = templates(:,:,:,[1,end]);
% templateG = m.templateG;

% ntemplates = size(templateG,4); %size(m,'templates'); ntemplates = ntemplates(4);
% save([savepath,filename(1:end-4),'reg.mat'],'R','Y','templates','-v7.3'); %,'shiftStruct'
save([savepath,filename],'R','Y','templates','-v7.3'); %,'shiftStruct'
mfile = matfile([savepath,filename],'Writable',true);

try
    movefile([filepath,filename],[filepath,'done/',filename])
catch
    mkdir([filepath,'done/'])
    movefile([filepath,filename],[filepath,'done/',filename])
end
disp(filename)


landmarkThreshold = 100;
template = templates(:,:,:,end).*(templates(:,:,:,end)>landmarkThreshold) + mask;



% finest correction (1)
options = NoRMCorreSetParms('d1',sizY(1),'d2',sizY(2),'d3',sizY(3)','grid_size',[10,20,10],'bin_width',10,'mot_uf',2,'us_fac',10,...
    'method',{'median','mean'},'overlap_pre',5,'overlap_post',5,'max_dev',[5,5,5],...
    'mem_batch_size',10,'use_parallel',false,...
    'init_batch',10,'iter',5,'max_shift',5,'shifts_method','cubic','upd_template',false);

Rth = R.*(R>landmarkThreshold);
clear R
Rth = Rth+repmat(mask,1,1,1,sizY(4));
[Rth,shifts] = normcorre_batch_even(Rth,options,template);
clear Rth

Y = apply_shifts(mfile.Y,shifts,options);
mfile.Y = Y;
mfile.options = options;
%clear Y

R = apply_shifts(mfile.R,shifts,options);
mfile.R = R;
%clear R


% finest correction (2)
options = NoRMCorreSetParms('d1',sizY(1),'d2',sizY(2),'d3',sizY(3)','grid_size',[7,15,7],'bin_width',10,'mot_uf',2,'us_fac',10,...
    'method',{'median','mean'},'overlap_pre',3,'overlap_post',3,'max_dev',[3,3,3],...
    'mem_batch_size',10,'use_parallel',false,...
    'init_batch',10,'iter',5,'max_shift',5,'shifts_method','cubic','upd_template',false);

Rth = R.*(R>landmarkThreshold);
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
clear R

