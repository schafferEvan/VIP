for i = 1:9
sourceFile = ['C:\Users\Axel-SCAPE\Documents\calibrations_10_6_17\SCAN\tiff_stacks\SCAN' num2str(i) '\unCorrected_SCAN' num2str(i) '.tiff'];
destinationFile = ['C:\Users\Axel-SCAPE\Documents\calibrations_10_6_17\SCAN\tiff_stacks\unCorrected_SCAN' num2str(i) '.tiff'];

movefile(sourceFile, destinationFile)
end
%%
clc 
close all
fileName = ['MAX_ALL_SCANS.tif'];
for i =1:9
   scan(:,:, i) = (imread(fileName, i));
end

%%

HR_timecourse = squeeze(sum(moviesub(:, :,:), 2));
final = cat(2, HR_timecourse, timecourses);
C = corrcoef(double(final));
vals = C(size(HR_timecourse, 2)+1:end, 1:size(HR_timecourse, 2));

for i = 1:size(vals, 1)
   peaks(i) =  find(vals(i, :) == max(vals(i, :)))
end
differences = diff(peaks)
keyboard
differences(end-1:end) = [];
differences = abs(differences);
stepSize_angular = info.daq.scanAngle/(info.daq.pixelsPerLine-2); % Remove 1 for flyback and another for the zero position. For a 100 um scan, with 102 pixels per line, there are 100 steps or 1 um
stepSize_um = 25; %um

cal_um_per_angle_mean = mean(stepSize_um./(stepSize_angular*differences));
cal_um_per_angle_std = std(stepSize_um./(stepSize_angular*differences));

%%
figure(1), colormap jet;
for i = 1:9
    imagesc(squeeze(scan(:, :, i)), [100 10000]);
    title(num2str(i));
    [x, y] = ginput(2);
    x = round(x);
    y = round(y);
    tmp = squeeze(scan(min(y):max(y), min(x):max(x), i));
    [a, b] = ind2sub(size(tmp), find(tmp == max(max(tmp))));
    a = round(mean(a));
    b = round(mean(b));
    xInd(i) = min(x)+b-1;
    yInd(i) = min(y)+a-1;
    
    hold on
    plot(xInd(i), yInd(i), 'ko')
    hold off
    pause(1)
end
differences = diff(yInd);
close all;
stepSize_um = 25;
cal_um_per_angle_mean = mean(stepSize_um./(stepSize_angular*differences))
cal_um_per_angle_std = std(stepSize_um./(stepSize_angular*differences))
