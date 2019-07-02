
% Clears everything
clear; clc; %close all;


% Select the directory in which the runs (scans) are stored
handles.info.directory = [uigetdir('D:\Data\','Select Parent Experiment directory to READ') '\'];
if 0 == handles.info.directory
    errordlg('You did not select a directory. Please choose again.','yeah')
    return;
end
handles.info.directory = handles.info.directory(1:end-1);
dataDirectory = handles.info.directory;
filesInDataDirectory = dir(handles.info.directory);

% Parse through all of the files in your data directory and compile a list
% of run names.
nameCounter=0;
for i = 1:length(filesInDataDirectory)
    if filesInDataDirectory(i).isdir ==1 & strcmp(filesInDataDirectory(i).name,'movies')==0;
        nameCounter=nameCounter+1;
        runNames{nameCounter} = filesInDataDirectory(i).name;
    end
end

% Open a list selection dialog box from which the user can load in a
% desired run
[SELECTION,OK] = listdlg('ListString',runNames)




% Start loading in that run
for qqq = SELECTION;
    scanName = runNames{qqq};
    infoFilePath = fullfile(dataDirectory, [scanName '_info.mat']);
    cd (dataDirectory)
    if exist(infoFilePath)
        % load info file
        disp('LOADING INFO FILE')
        load(infoFilePath);
        disp(info.daq)
        
        % Scake data set size by spatial bin factor
        binFactor = info.camera.binFactor;          % Camera Bin Number
        numDepths = 250;%info.camera.yROI / binFactor;   % Number of Depth Pixels
        numLatPix = 800;%info.camera.xROI / binFactor;     % Number of Lateral Pixels
        %         numDepths = info.camera.yROI;
        %         numLatPix = info.camera.xROI;
        scanRate = info.daq.scanRate;               % Volumetric Scan Rate (VPS)
        
        
        % Given in microns/pixel for depth, lateral and scan dimensions (in that
        % order)
        % Set the conversion factors based on objective used
%         presentDirectory = pwd;
%         cd C:\Users\Andor\Documents\Matt\LSIPT_galvo_drive_GUI;
%         load savedConversionFactors.mat
%         cd(presentDirectory)
%         objective = info.objective;
%         
%         
%         if (isempty(strfind(objective, '10x')) == 0)
%             cFact = savedConversionFactors(1);
%             conversionFactors(1) = cFact.depth_umPerPix;     % Does change as a function of objective magnification
%             conversionFactors(2) = cFact.lateral_umPerPix;   % Does change as a function of objective magnification
%             conversionFactors(3) = 1/(cFact.scan_VoltsPerUm*info.daq.pixelsPerLine/info.daq.scanAngle); % Doesn't change as a function of objective magnification (Volts to galvo/um)alvo/um)
%         elseif(isempty(strfind(objective, '20x')) == 0)
%             cFact = savedConversionFactors(2);
%             conversionFactors(1) = cFact.depth_umPerPix;     % Does change as a function of objective magnification
%             conversionFactors(2) = cFact.lateral_umPerPix;   % Does change as a function of objective magnification
%             conversionFactors(3) = 1/(cFact.scan_VoltsPerUm*info.daq.pixelsPerLine/info.daq.scanAngle); % Doesn't change as a function of objective magnification (Volts to galvo/um)
%         end
%         
        
        %conversionFactors = 1./[0.4724/binFactor .5512/binFactor 0.0028*info.daq.pixelsPerLine/info.daq.scanAngle];
        %conversionFactors = 1./[1.414/binFactor 1.0709/binFactor 0.0031*info.daq.pixelsPerLine/info.daq.scanAngle];
        
        disp('LOADING METADATA')
        % Creates file path of zyla metadata file (.ini)
        zylaInfoFilePath = fullfile(dataDirectory, scanName, 'acquisitionmetadata.ini');
        FID = fopen(zylaInfoFilePath, 'r');
        zylaMetaData = fread(FID, '*char')';
        fclose(FID);
        % Reads the number of frames in a spool file
        startIndex = strfind(zylaMetaData, 'ImagesPerFile') + length('ImagesPerFile = ');
        numFramesPerSpool = str2double(zylaMetaData(startIndex:end));
        % Reads the image size in the spool file
        startIndex = strfind(zylaMetaData, 'ImageSizeBytes') + length('ImageSizeBytes = ');
        ImageSize = str2double(zylaMetaData(startIndex : startIndex + 7));
        % Reads in AOIStride
        startIndex = strfind(zylaMetaData, 'AOIStride') + length('AOIStride = ');
        AOIStride = str2double(zylaMetaData(startIndex : startIndex+3))/2;
        
        
        numFrames = length(info.daq.scanWaveform);  % Total number of frames acquired
        numScans = info.daq.numberOfScans;          % Total number of volumes acquired
        framesPerScan = floor(info.camera.framerate/info.daq.scanRate);        % Number of frames per volume
        

        % add extra columns for 2 buffer rows (why?)
        if (2 == info.camera.binFactor)
            numColumns = numDepths + 1;
            numRows = ImageSize / 2 / numColumns;
            if(mod(numRows, 1) ~= 0)
                numColumns = numDepths+2;
                numRows = ImageSize/2/numColumns;
            end
        else
            numColumns = numDepths + 2;
            numRows = AOIStride;
        end
        
        numFramesAcquired = round(str2num(info.camera.kineticSeriesLength));
        
        % Creates a list of spool file names that exist for a given run
        disp('CREATING SPOOL FILE NAMES')
        spoolCounter=0;
        for i = 1:floor(numFramesAcquired/(numFramesPerSpool))
            spoolCounter=spoolCounter+1;
            temp = i;
            for j = 1:10
                a(i, j) = mod(temp, 10^j)/(10^(j-1));
                temp = temp-mod(temp, 10^j);
            end
            tempName = mat2str(a(i, :));
            tempName = tempName(2:end-1);
            tempName = tempName(find(tempName ~= ' '));
            tempName = [tempName 'spool.dat'];
            namesOut{spoolCounter} = tempName;
        end
        namesOut = [{'0000000000spool.dat'} namesOut];
        filesToLoad = namesOut
        clear namesOut tempName temp a spoolCounter
        
        cd(dataDirectory);
        cd(scanName);
        disp('ALLOCATING MEMORY')
        tic
        moviesub = uint16(0);
        moviesub(numRows, numColumns, numFramesAcquired) = 0;
        toc
        spoolFileSize = numRows*numColumns*numFramesPerSpool;
        exitCondition = 0;
        fileCounter = 1;
        h = waitbar(0, 'Please Wait');
        movieSize = prod(size(moviesub));
        
        while exitCondition == 0
            
            firstIndex = (fileCounter-1)*spoolFileSize +1;
            lastIndex = (fileCounter*spoolFileSize);
            if (firstIndex > movieSize) % All of these pixels will be dark images
                exitCondition = 1;
            elseif (lastIndex > movieSize) % Some of these pixels will be dark images
                dataLeftToFill = [firstIndex:movieSize];
                if (length(dataLeftToFill) ~= 0)
                    fileToRead = filesToLoad{fileCounter};
                    fid = fopen(fileToRead, 'r');
                    X = fread(fid, 'uint16=>uint16')';
                    fclose(fid);
                    moviesub(dataLeftToFill) = ...
                        X(1:length(dataLeftToFill));
                end
                exitCondition = 1;
            else % All of these pixels will be moviesub
                fileToRead = filesToLoad{fileCounter};
                fid = fopen(fileToRead, 'r');
                moviesub(firstIndex:lastIndex) = ...
                    fread(fid, 'uint16=>uint16')';
                fclose(fid);
            end
            waitbar(lastIndex/movieSize)
            fileCounter = fileCounter +1;        
        end
        
        
        darkImage = moviesub(:, :, numFrames+1:end);
        meanDark = squeeze(mean(mean(darkImage, 1), 2));
        toKeep = find(meanDark < 110);
        darkImage = darkImage(:, :, toKeep);
        moviesub(:, :, numFrames+1:end) = [];
        
        close (h)
    else
        keyboard
    end
    clearvars -except qqq runNames SELECTION darkImage binFactor moviesub dsf info scanName dataDirectory
    %FOV = size(moviesub).*[conversionFactors(2) conversionFactors(1) conversionFactors(3)];
    moviesub(:, end-1:end, :, :) = [];
    darkImage(:, end-1:end, :, :) = [];

    figure(1)
    latRange = [400 500];
    A = squeeze(max(moviesub(latRange(1):latRange(2), :, :), [], 3));
    imagesc(A)
    pause(1)
    index = num2str(find(qqq==SELECTION));
    disp([info.scanName ': ' num2str(index)])
    pause(1)
    cd(dataDirectory)
    savefile = ['scanMIP_full' index '.mat']
    save(savefile, 'A');
    
end

cd(dataDirectory)
disp('FINISHED SCANS IN DATA');
%%
numFiles = length(SELECTION);
dataSet = zeros(size(A, 1), size(A, 2), numFiles);
clear A
for i = 1:length(SELECTION)
    loadFile = ['scanMIP_full' num2str(i) '.mat'];
    load(loadFile);
    dataSet(:, :, i) = A;
    clear A;
end

dataSet = squeeze(max(dataSet, [], 3));
dataSet = dataSet';

%%
figure(1), colormap gray
imagesc(dataSet);
xlabel('Lateral');
ylabel('Depth');

for i = 1:length(SELECTION)
   [x, y] = ginput(2)
   title('Bound the bead you care about')
   
   x = round(x), y = round(y);
   bead = dataSet(min(y):max(y), min(x):max(x));
   
   profile = squeeze(max(bead, [], 2));
   index = find(profile == max(profile))
   depths(i) = index+min(y);
end

depths
differences = diff(depths)
keyboard

stepSize = 25.4;
depth_um_per_pix_mean = mean(stepSize./differences)
depth_um_per_pix_std = std(stepSize./differences)
