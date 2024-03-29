function make_rawScaled_movies_twoColor(codePath, experimentFolder, savePath)
% skips frames to quickly make movie containing sample frames from entire
% experiment (primarily to check motion correction)

addpath(genpath(codePath))
infoFile = dir([experimentFolder,'info/f*.mat']);
%pause(60*150)

% if ismac
%     addpath(genpath('~/Dropbox/_AxelLab/matlab/calcium-signal-extraction/'))
%     addpath(genpath( '~/Dropbox/_code/pipeline/'))
%     experimentFolder = '/Volumes/dataFast2/2019_05_14_57C10_NLS6f_CrzDh44_red_Su_walk/';%'/Volumes/SCAPEdata1/finalData/2019_03_12_Nsyb_NLS6s_Su/fly2/';%'/Volumes/dataFast/habaRegistered/2019_02_14_Nsyb_NLS6s_Su/';
%     infoFile = dir([experimentFolder,'info/*.mat']);%'/Volumes/dataFast/habaRegistered/2018_08_24_odor/mats/new/fly3run2/'; %'/Users/evan/Desktop/hungerRaw/';%'/Volumes/SCAPEdata1/scratchData/2018_08_01_IRtest/matfiles/registered/';%'/Volumes/data/_scape/data/_outMats/'; %
%     if length(infoFile)>1; infoFile=infoFile(end); end
%     savePath = '~/Dropbox/_AxelLab/_data/_scape/movies/';
% elseif isunix
%     addpath(genpath('/home/analysis-pc/00_Analysis/motion_correction/'))
%     addpath(genpath('/home/analysis-pc/00_Analysis/calcium-signal-extraction/'))
%     experimentFolder = '/home/analysis-pc/rawData/20171025_nSyb_fly2_reg/registered/reregistered/';
%     savePath = '/home/analysis-pc/rawData/movies/';
% else
%     addpath(genpath('E:\Dropbox\GitHub\calcium-signal-extraction'))
%     addpath(genpath('E:\Dropbox\GitHub\eftyMotionCorr'))
%     experimentFolder = 'D:\SCAPEdataworkingfolder\_outMats\';
%     savePath = 'D:\SCAPEdataworkingfolder\movies\';
% end

% if isfolder([experimentFolder,'Yproj/']); load([experimentFolder,'Yproj/YprojFit.mat'],'YprojFit');scaling=YprojFit/max(YprojFit);
% else; scaling = [];
% end
scaling = [];
[trials,~, ~, runNum] = sortExperimentDirectory(experimentFolder,'reg');

params.trialName = [trials(end).name(1:end-4),'all.avi']; %[trials(end).name(1:runLoc-1),'all.avi'];
params.savePath = savePath;
params.concatenate = true;
load([experimentFolder,'info/',infoFile(1).name]);
params.Ttot = 0;
params.acqRate = round(info.daq.scanRate); %10; % volumes per second
chunkSize = 101;

% select maxY to set color scale
m = matfile([experimentFolder,trials(1).name]);
try
    sz = size(m,'Ysmall');
    Y = single( m.Ysmall(:,:,:,sz(4)-50:sz(4) ) );
catch
    sz = size(m,'Y');
    Y = single( m.Y(:,:,:,sz(4)-50:sz(4) ) );
end
params.maxY = max(Y(:));

try
    params.xum = sz(1)*info.GUIcalFactors.x_umPerPix;
catch
    %warning('no x_um field available')
    params.xum = sz(1)*info.GUIcalFactors.xK_umPerVolt*info.daq.scanAngle/(info.daq.pixelsPerLine-1);
end
params.yum = sz(2)*info.GUIcalFactors.y_umPerPix;
params.zum = sz(3)*info.GUIcalFactors.z_umPerPix;


getPID = false; %true;
if getPID
    fid = fopen('/Volumes/SCAPEdata1/rawData/2018_08_24_NsybNLS_odors/fly3_run2_4reps_seq_stim_data.bin','r');
    stimData = fread(fid,[8,inf],'double');
    fclose(fid);
    pidtime = squeeze(stimData(1,:));
    pidval = stimData(6,:);% 6 for PID, 8 for odor sensor
    onTh = 0.125;
    onMax = max(pidval);
    pidResampled = interp1(pidtime, pidval, 30+1/params.acqRate:1/params.acqRate:pidtime(end));
    pidResampled = (pidResampled-onTh).*(pidResampled>onTh)/(onMax-onTh);
end

rTh = 250;
mask = m.R(:,:,:,1)>rTh;

disp(['making BIG movie of ',num2str(length(trials)),' trials'])
params.kTot=0;
%% load in registered data and find ROIs
for i=1:length(trials)
    if i>1; if runNum(i)~=runNum(i-1); params.Ttot = 0; end; end
    
    trialPath = [experimentFolder,trials(i).name];
    display(['file ',num2str(i),': ',trials(i).name])
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     x = h5info(trialPath,'/mov');
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     %sz = x.Dataspace.Size;
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     Y = bigread2(trialPath);
    m = matfile(trialPath);
    try
        sz = size(m,'Ysmall');
    catch
        sz = size(m,'Y');
    end
    pts = [1:chunkSize:sz(end),sz(end)+1]; 
    %if pts(end)<sz(end); pts = [pts,sz(end)]; end
    for k=1:length(pts)-1
        try
            Y = single( m.Ysmall(:,:,:,pts(k):pts(k+1)-1) );
        catch
            Y = single( m.Y(:,:,:,pts(k):pts(k+1)-1) );
        end
        szy = size(Y);
        if ~isempty(scaling)
            dvec = scaling(params.kTot+(1:szy(4)));
            for t=1:szy(4)
                Y(:,:,:,t)=Y(:,:,:,t)/(dvec(t)*params.maxY);
            end
        else
            Y = Y/params.maxY;
        end
        
%         % manual extra cropping to match caiman movie
%         Y = Y(8:102,6:254,5:114,:);
%         Y = Y(1:91,1:246,:,:);
%         
%         load(['/Users/evan/Dropbox/temp/_0110temp/','alignedBehavAndStim.mat'],'alignedStim');
%         stimOn = alignedStim;
%         stimOn = [stimOn,ones(1,4)];
%         %stimOn = [2600 3930]; %pidResampled( round(params.Ttot*params.acqRate)+(1:szy(4)) ); %[inf,inf]; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%[300,900];
        %         stimOn = stimOn-params.Ttot;
        %         if stimOn(2)>=1
        %             if stimOn(1)<1; stimOn(1)=1; end
        %         else
        stimOn = [inf,inf];
        %         end

        params.vid = plot4Dproj_ES(Y, mask, [szy(1),szy(2),szy(3)],params,stimOn);
        params.Ttot = params.Ttot + .1*round(10*szy(4)/params.acqRate);
        params.kTot=params.kTot+szy(4);
        close(gcf)
    end
end
close(params.vid)
end
