clc, clear, close all

for i = 2:11
   fileName = ['depth' num2str(i-1) '.tif'];
   lateral(:,:, i+1) = imread(fileName);
end

figure(1)
timecourses = squeeze(sum(lateral(:, :, :), 1));
plot(timecourses(:, 1:end));
%%
close all
for i = 1:size(timecourses, 2)
    plot(timecourses(:, i))
    [x, y] = ginput(2);
    x = round(x);
    peakIndex(i) = min(x)+round(mean(find(timecourses(min(x):max(x), i)...
        == max(timecourses(min(x):max(x), i)))))
end
differences = diff(peakIndex);
differences(1) = [];
stepSize = 25.4; %um
lateralCal_mean = mean(stepSize./differences)
lateralCal_std = std(stepSize./differences)
