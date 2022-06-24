function D = ImaGIN_BadChannel(S)
% Detect bad channels indices with a trained classifier.
%
% -=============================================================================
% This function is part of the ImaGIN software: 
% https://f-tract.eu/
%
% This software is distributed under the terms of the GNU General Public License
% as published by the Free Software Foundation. Further details on the GPLv3
% license can be found at http://www.gnu.org/copyleft/gpl.html.
%
% FOR RESEARCH PURPOSES ONLY. THE SOFTWARE IS PROVIDED "AS IS," AND THE AUTHORS
% DO NOT ASSUME ANY LIABILITY OR RESPONSIBILITY FOR ITS USE IN ANY CONTEXT.
%
% Copyright (c) 2000-2018 Inserm U1216
% =============================================================================-
%
% Authors: Viateur Tuyisenge & Olivier David

% Name of input iEEG file .dat/.mat (cropped)
if (nargin >= 1) && isfield(S, 'dataset') && ~isempty(S.dataset)
    FileIn = S.dataset; 
else
    FileIn = spm_select(1, '\.mat$', 'Select EEG mat file');
    if isempty(FileIn)
        return;
    end
end

% Directory where iEEG file .dat/.mat with badchannels indices will be stored
if (nargin >= 1) && isfield(S, 'FileOut') && ~isempty(S.FileOut)
    [badDir, FileOut, ~] = fileparts(S.FileOut);
    isNewFile = 1;
else
    [badDir, FileOut, ~] = fileparts(FileIn);
    isNewFile = 0;
end

% Directory where training set is stored (if not defined, use the toolbox folder)
if (nargin >= 1) && isfield(S, 'trainBase') && ~isempty(S.trainBase)
    trainDir = S.trainBase;
else
    trainDir = fileparts(mfilename('fullpath'));
end

if (nargin >= 1) && isfield(S, 'toCheck') && ~isempty(S.toCheck)
    toCheck = S.toCheck;
else
    toCheck = 'toCheck';
end

% Detect bad channels indices and save summary images
% Load the cropped meeg object
D = spm_eeg_load(FileIn);

chanLbs = D.chanlabels;
elec = sensors(D,'eeg');  
pos  = elec.elecpos; 
Sens = elec.label;
idxNaN = find(isnan(pos(:,1)));
% if badchanel file exists, "text file ending with _bIdx"
% no computation is needed, load and put them in spm object.
if exist(fullfile(badDir, [FileOut,'_bIdxChecked.txt']),'file') == 2 && strcmpi(toCheck,'OK')    
    bIdxAuto    = load(fullfile(badDir, [FileOut,'_bIdx.txt']));   % txt file of indices computed by machine learning
    bIdxChecked = load(fullfile(badDir, [FileOut,'_bIdxChecked.txt'])); % txt file of indices manually written during quality control
    if ~isempty(bIdxChecked)
        bIdx = [bIdxChecked(:);idxNaN(:)];
    else
        bIdx = [bIdxAuto(:);idxNaN(:)];
    end
    bIdx     = sort(unique(bIdx));
    interIdx = intersect(bIdx,idxNaN);
    NaNbIdx  = setdiff(bIdx,idxNaN);
    ixNaN    = setdiff(idxNaN, bIdx);
    fprintf('MESSAGE: used bIdx previously validated!')
else
    % Do all computations    
    clear S2;
    S2.FileName = FileIn;
    if ~isempty(D.events) && ~isempty(D.events.type) && ismember('Stim', {D.events.type})
        S2.InterpolationFilter = 1;
    else
        S2.InterpolationFilter = 0;
    end
    
    T = ImaGIN_FeatureSEEG(S2);
    
    % If the trained classifier is not available: compute it
    trainedFile = fullfile(trainDir, 'ImaGIN_trainedClassifier.mat');
    if ~exist(trainedFile, 'file')
        % Get train base file
        trainBaseFile = fullfile(trainDir, 'ImaGIN_trainBaseFeatures.mat');
        if ~exist(trainBaseFile, 'file')
            error(['File ' trainBaseFile ' was not found.']);
        end
        % Load the train base
        trainBase = load(trainBaseFile);
        % Train the classifier
        trainedClassifier = ImaGIN_trainClassifier(trainBase.predictors, trainBase.response);
        % Save results
        save(trainedFile, '-struct', 'trainedClassifier');
        % Otherwise, load the trained classifier
    else
        trainedClassifier = load(trainedFile);
    end
    
    % Predict new dataset
    channelClass = trainedClassifier.predictFcn(T(:,2:8));
    % Get list of detected bad channels
    bIdx = T.noIdx(strcmp(channelClass, 'Bad'));
    
    %%
    % In case disconnected electrode doesn't have stimulation artefact
    % specific for some FTRACT datasets
    
   
    
    crFname = D.fname;
    crFname = strrep(crFname,'welectrodes_','');
    undsc   = strfind(crFname,'_'); % Ap12_3mA_1Hz_1050us becomes Ap01-Ap02_3mA_1Hz_1050us with new convention
    idxDush = strfind(crFname,'-');
    Pulse   = strfind(crFname,'us');
    Amp     = strfind(crFname,'mA');
    Frq     = strfind(crFname,'Hz');
    if numel(undsc) == 4 && ~isempty(Pulse)&& ~isempty(Amp) && ~isempty(Frq)
        if ~isempty(idxDush)
            ch1 = crFname(1:idxDush(1)-1);
            ch2 = crFname(idxDush(1)+1:undsc(1)-1);
        else
            ch1 = '';
            ch2 = '';
        end
        
        iLastLetter1 = find(~ismember(ch1, '0123456789'), 1, 'last');
        iLastLetter2 = find(~ismember(ch2, '0123456789'), 1, 'last');
        % If there is not electrode label or no index: wrong naming
        if isempty(iLastLetter1) || (iLastLetter1 == length(ch1))
            return;
        end
        
        if isempty(iLastLetter2) || (iLastLetter2 == length(ch2))
            return;
        end
        % Find channel index
        chLabel1 = ch1(1:iLastLetter1);
        chInd1 = str2num(ch1(iLastLetter1+1:end));
        chLabel2 = ch2(1:iLastLetter2);
        chInd2 = str2num(ch2(iLastLetter2+1:end));
        
        ch1 = sprintf('%s%0d', chLabel1, chInd1);
        ch2 = sprintf('%s%0d', chLabel2, chInd2);
        
        ch1_2d = sprintf('%s%02d', chLabel1, chInd1);
        ch2_2d = sprintf('%s%02d', chLabel2, chInd2);
        
        chlb = {ch1,ch2};
        
        chlb_2d = {ch1_2d,ch2_2d};
        if ~isempty(chlb) || ~isempty(chlb_2d)
            chInd1  = find(strcmp(chanLbs,ch1));
            chInd2  = find(strcmp(chanLbs,ch2));
            if isempty(chInd1)
                chInd1  = find(strcmp(chanLbs,ch1_2d));
            end
            if isempty(chInd2)
                chInd2  = find(strcmp(chanLbs,ch2_2d));
            end
            chInd = [chInd1(:);chInd2(:)];
            
            if ~isempty(chInd)
                if isempty(find(any(bIdx==chInd(1)), 1))
                    bIdx(end+1) = [chInd(1)];
                end
                if numel(chInd) > 1
                    if isempty(find(any(bIdx==chInd(2)), 1))
                        bIdx(end+1) = chInd(2);
                    end
                end
                bIdx = sort(bIdx);
            end
        end
    end
    bIdx = bIdx(:);    
    interIdx = intersect(bIdx,idxNaN);
    NaNbIdx  = setdiff(bIdx,idxNaN);
    ixNaN    = setdiff(idxNaN,bIdx);
    if ~isempty(idxNaN)
        bIdx = [bIdx(:);idxNaN(:)];
        bIdx = sort(unique(bIdx));
    end
    
    %%
    % Save bad channel indices in .txt file
    
    badIdxFile = fopen(fullfile(badDir, [FileOut, '_bIdx.txt']), 'w');
    badchaFile = fopen(fullfile(badDir, [FileOut, '_bChans.txt']), 'w');
    for i = 1:length(bIdx)
        fprintf(badchaFile, '%d %s\n', bIdx(i), char(chanLbs{bIdx(i)}));
        fprintf(badIdxFile, '%d \n', bIdx(i));
    end
    fclose(badchaFile);
    fclose(badIdxFile);
    
    Lia = ismember(T.noIdx,bIdx);
    channelClass(Lia) = {'Bad'};
    Tnew = [T channelClass];
    Tnew.Properties.VariableNames{'Var9'} = 'Note';
    csvfilename = fullfile(badDir, [FileOut, '.csv']); % Save feature table & badchannels indices
    writetable(Tnew,csvfilename,'Delimiter',',');
    
end

try
    monoRecordings = fopen(fullfile(badDir, ['recordings_monopolar_', FileOut, '.txt']), 'w'); % export monopolar recording channels
    for i = 1:length(Sens)
        if isempty(find(strcmp(Sens{i}, chanLbs(bIdx)),1))
            if isempty(regexpi(char(Sens{i}),'ecg'))
                fprintf(monoRecordings, '%s\n', char(Sens{i}));
            end
        end
    end
    fclose(monoRecordings);
catch exception 
    throw(exception)
end

% Add badchannel index in meeg object
if ~isempty(bIdx)
    D = badchannels(D,bIdx,1);
end
% Create new file in output (otherwise, simply update the input file)
if isNewFile
    Dbad = clone(D, fullfile(badDir, FileOut), [D.nchannels D.nsamples D.ntrials]); % save meeg with badchannel indices in Badchannel directory
    Dbad(:,:,:) = D(:,:,:);
    % Save modified list of bad channels
    save(Dbad);
else
    % Save modified list of bad channels
    save(D);
end

% Channel plots and ScreenShots
figDir = fullfile(badDir, 'ScreenShot');
if ~exist(figDir, 'dir')
    mkdir(figDir);
end

% close all;
Size = 8;  % Number of channels per screenshot
n_c  = size(D,1);
if n_c >= Size
    tmp = floor(n_c/Size);
    for i2 = 1:tmp
        figure(i2);
        set(gcf,'Position',[629 -17 702 1101])
        for i3 = 1:Size  
            if ~ismember(chanlabels(D,i3+(i2-1)*Size),D.csv_labels)
                color = 'g'; % Channel was not matched with the .csv
            elseif intersect(i3+(i2-1)*Size,ixNaN) == i3+(i2-1)*Size
                color = 'b'; % Channel was matched, is good, has NaN coordinates            
            elseif intersect(i3+(i2-1)*Size,interIdx) == i3+(i2-1)*Size 
                color = 'm'; % Channel was matched, is bad, has NaN coordinates 
            elseif intersect(i3+(i2-1)*Size,NaNbIdx) == i3+(i2-1)*Size  
                color = 'r'; % Channel was matched, is bad, has coordinates
            else
                color = 'k'; % Channel was matched, is good, has coordinates
            end
            subplot(Size,1,i3)
            plot(time(D),D(i3+(i2-1)*Size,:),color);
            ylabel([num2str(i3+(i2-1)*Size) ' : ' D.chanlabels{i3+(i2-1)*Size}])
            if i3 == 1
                figName = char(strcat(FileOut,'_',num2str(i3+(i2-1)*Size),'-', ...
                    num2str(i2*Size)));
                title(figName,'interpreter','none');
            end
            axis tight
        end
        zoom on
        fig = figure(i2);
        print(fig,fullfile(figDir,figName),'-dpng'); %ScreenShot
        close;
    end
    
    rmd = size(D,1) - tmp*Size;
    if rmd ~= 0
        figure(tmp + 1)
        set(gcf,'Position',[629 -17 702 1101])
        for i4 = 1:rmd             
            if ~ismember(chanlabels(D,i3+(i2-1)*Size+i4),D.csv_labels)
                color = 'g'; % Channel was not matched with the .csv
            elseif intersect(i3+(i2-1)*Size+i4,ixNaN) == i3+(i2-1)*Size+i4
                color = 'b'; % Channel was matched, is good, has NaN coordinates            
            elseif intersect(i3+(i2-1)*Size+i4,interIdx) == i3+(i2-1)*Size+i4 
                color = 'm'; % Channel was matched, is bad, has NaN coordinates 
            elseif intersect(i3+(i2-1)*Size+i4,NaNbIdx) == i3+(i2-1)*Size+i4  
                color = 'r'; % Channel was matched, is bad, has coordinates
            else
                color = 'k'; % Channel was matched, is good, has coordinates
            end
            subplot(rmd,1,i4)
            plot(time(D),D(i3+(i2-1)*Size+i4,:),color);
            ylabel([num2str(i3+(i2-1)*Size+i4) ' : ' D.chanlabels{i3+(i2-1)*Size+i4}])
            if i4 == 1
                figName = char(strcat(FileOut,'_',num2str(i3+(i2-1)*Size + 1),'-',num2str(n_c)));
                title(figName,'interpreter','none');
            end
            axis tight
        end
        zoom on
        fig = figure(i2+1);
        print(fig, fullfile(figDir,figName), '-dpng');
        close
    end
else
    figure(1)
    set(gcf,'Position',[629 -17 702 1101])
    for i5 = 1:n_c
        if ~ismember(chanlabels(D,i5),D.csv_labels)
            color = 'g'; % Channel was not matched with the .csv
        elseif intersect(i5,ixNaN) == i5
            color = 'b'; % Channel was matched, is good, has NaN coordinates 
        elseif intersect(i5,interIdx) == i5
            color = 'm'; % Channel was matched, is bad, has NaN coordinates
        elseif intersect(i5,NaNbIdx) == i5
            color = 'r'; % Channel was matched, is bad, has coordinates
        else
            color = 'k'; % Channel was matched, is good, has coordinates
        end         
        subplot(n_c,1,i5)
        plot(time(D),D(i5,:),color);
        ylabel([num2str(i5) ' : ' D.chanlabels{i5}])
        if i5 == 1
            figName = char(strcat(FileOut,'_',num2str(1),'-',num2str(n_c)));
            title(figName,'interpreter','none');
        end
        axis tight
    end
    zoom on
    fig = figure(1);
    print(fig, fullfile(figDir,figName), '-dpng');
    close
end
set_final_status('OK')
    
