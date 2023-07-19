

function compile_motion_metrics(codePath, experimentFolder, isRegistered, color_quants, z_th)
% color_quant = 0.95
% z_th = 100

addpath(genpath(codePath))

if isRegistered
    [trials, trialOrder] = sortExperimentDirectory(experimentFolder,'reg',false);
else
    [trials, trialOrder] = sortExperimentDirectory(experimentFolder,'',false);
end

% more extensive check of whether files are corrupt
triggerError = false;
for ii=1:length(trials)    
    trialPath = [experimentFolder,trials(trialOrder(ii)).name];
    lastwarn('') % Clear last warning message    
    try
        m = matfile(trialPath);
        sy = size(m,'Y');
        sr = size(m,'R');
        m.Y(sy(1),sy(2),sy(3),sy(4));
        m.R(sr(1),sr(2),sr(3),sr(4));
        m.Y(sy(1),sy(2),sy(3),1);
        m.R(sr(1),sr(2),sr(3),1);
    catch
        triggerError = true;
        warning([trials(trialOrder(ii)).name,' is corrupted']);
    end
    [~, warnId] = lastwarn;
    if strcmp(warnId,'MATLAB:whos:UnableToRead')
        disp([trials(trialOrder(ii)).name,' is PARTIALLY corrupted']);
    end
end
if triggerError
    error('Aborting due to corrupted files')
end
if nargin<4
    color_quants = [.75,.999];
end
if nargin<5
    disp('z_th unused')
    %z_th = sr(3);
end

%% sum image calc
cR=[];
for ii=1:length(trials)
    
    trialPath = [experimentFolder,trials(trialOrder(ii)).name];
    display(['file ',num2str(ii),': ',trials(trialOrder(ii)).name])
    m = matfile(trialPath);
    if isRegistered
        R = m.R;
    else
        R = m.Y;
    end
    
    % ref
    q = quantile(R(:),color_quants);
    ref_img = mean(R,4);
    ref_img(ref_img>q(2))=q(2);
    ref_img = ref_img-q(1);
    ref_img(ref_img<0)=0;

    % sample
    q = quantile(R(:),color_quants);
    %R = R(:,:,1:z_th,:);
    R(R>q(2))=q(2);
    R = R-q(1);
    R(R<0)=0;
    
    
    [cRtmp,~,~] = motion_metrics_on_static_ref(R,ref_img);
    cR = [cR;cRtmp];
     
end

mkdir([experimentFolder,'Yproj'])
if isRegistered
    save([experimentFolder,'Yproj/motion_metrics.mat'],'cR')
else
    save([experimentFolder,'Yproj/motion_metrics_raw.mat'],'cR')
end
