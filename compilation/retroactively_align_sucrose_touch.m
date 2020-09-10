
function retroactively_align_sucrose_touch(codePath, traceFolder)
% sucrose touch was not part of original behavior log. After generating
% timestamp for this event (in behavior frames), this function converts
% timestamp to units of seconds preceding first stim


if ~nargin
    addpath(genpath('..'))
    base_folder = '/Volumes/SCAPEdata1/finalData/';
    exp_id = '2019_04_18_Nsyb_NLS6s_Su_walk/';
    flynum = 'fly2';
    traceFolder = [base_folder, exp_id, flynum, '/Yproj/'];
    %behaviorDataDir = [];
else
    addpath(genpath(codePath))
    flynum = traceFolder(end-10:end-7);
end

% get beh cam fps
rawfile = [traceFolder, 'behavior/', flynum, '_run2.mat'];
% rawfile = [base_folder, exp_id, flynum, '/Yproj/behavior/', flynum, '_run2.mat'];
load( rawfile, 'traces' )
t1 = find(diff(traces.isStimOn)>0);
alignedfile = [traceFolder, 'alignedBehavAndStim.mat'];
%alignedfile = [base_folder, exp_id, flynum, '/Yproj/', 'alignedBehavAndStim.mat'];
load( alignedfile, 'alignedBehavior' )
t2 = find(diff(alignedBehavior.stim)>0);
beh_cam_fps = mean(diff(t1')./diff(alignedBehavior.timeTot(t2)));

% load touch time file
% touch_base = '/Volumes/SCAPEdata1/scapeBehavior/';
% touchfile = [touch_base, exp_id, flynum, '_run2_touch.mat'];
touchfile = [traceFolder, 'behavior/', flynum, '_run2_touch.mat'];
touch = load( touchfile );
touch_time = find(touch.traces.isStimOn>0, 1);

% relate touch time to first logged stim, in beh cam frames
beh_cam_touch_frame_dur = t1(1)-touch_time;

% compute touch timestamp in correctly aligned time
touch_timestamp = - beh_cam_touch_frame_dur/beh_cam_fps;

outfile = [traceFolder, 'sucrose_touch_aligned.mat'];
% outfile = [base_folder, exp_id, flynum, '/Yproj/', 'sucrose_touch_aligned.mat'];
save(outfile, 'touch_timestamp')

% %plot
% tc = zeros(size(alignedBehavior.stim));
% [~,v] = min(abs(alignedBehavior.timeTot-(round(alignedBehavior.timeTot(t2(1))+touch_timestamp)))
% tc(v:end)=1;
% figure; plot(alignedBehavior.timeTot, alignedBehavior.stim)
% hold all; plot(alignedBehavior.timeTot, tc)
