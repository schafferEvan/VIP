

function compile_motion_metrics(codePath, experimentFolder, isRegistered)

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

%% sum image calc
for ii=1:length(trials)
    
    trialPath = [experimentFolder,trials(trialOrder(ii)).name];
    display(['file ',num2str(ii),': ',trials(trialOrder(ii)).name])
    m = matfile(trialPath);
    if isRegistered
        R = m.R;
    else
        R = m.Y;
    end
    
    if ii==1
        [cR,ref_img,~] = motion_metrics_on_static_ref(R);
    else
        [cRtmp,~,~] = motion_metrics_on_static_ref(R,ref_img);
        cR = [cR,cRtmp];
    end
     
end

mkdir([experimentFolder,'Yproj'])
save([experimentFolder,'Yproj/motion_metrics.mat'],'cR')
