function D = ImaGIN_BadChannelQualityControl(S)
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

% FileIn: path linking to meeg object 
FileIn = S.dataset;
try
    D = spm_eeg_load(FileIn); % Load the cropped meeg object
catch
    FileIn = spm_select(1, '\.mat$', 'Select data file');
    D=spm_eeg_load(FileIn);
end
[badDir,cutName,~] = fileparts(FileIn);
bPrefix = strcat(badDir,'/',cutName); 

bIdxAuto    = load(strcat(bPrefix,'_bIdx.txt'));       % txt file of indices computed by machine learning
bIdxChecked = load(strcat(bPrefix,'_bIdxChecked.txt')); % txt file of indices manually written during quality control
if ~isempty(bIdxChecked)
    bIdx = bIdxChecked(:);
else
    bIdx = bIdxAuto(:);
end

elec= sensors(D,'eeg');  % add channels without positions into bad channels
pos = elec.elecpos; 
Sens = elec.label;
chanLbs = D.chanlabels;
idxNaN = []; % find(isnan(pos(:,1))); We now assume channels without coordinates are not bad (CRFs older than 2009 do not have T1pre coordinates)
if ~isempty(idxNaN)
    bIdx = [bIdx(:);idxNaN(:)];
    bIdx = sort(unique(bIdx));
end
T = readtable(strcat(bPrefix,'.csv'));
tNote = cell(size(T,1),1);
tNote(:) = {'Good'};
tNote(ismember(T.noIdx,bIdx)) = {'Bad'};
T(:,9) = tNote;
csvfilename = strcat(bPrefix,'.csv');
writetable(T,csvfilename,'Delimiter',',');

badchaFile = fopen(strcat(bPrefix,'_bChans.txt'), 'w');

if ~isempty(bIdx)
    for i = 1:length(bIdx)
        fprintf(badchaFile, '%d %s\n', bIdx(i), char(chanLbs{bIdx(i)}));
    end    
end
fclose(badchaFile);

try
    monoRecordings = fopen(fullfile(badDir, ['recordings_monopolar_', cutName, '.txt']), 'w'); % export monopolar recording channels
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

tmpIdx = D.badchannels; % old badchannels indices
D = badchannels(D,tmpIdx,0); % reset meeg object bad indices

if ~isempty(bIdx)
    D = badchannels(D, bIdx,1); % add new badchannel indeces in meeg object
end
Dbad = clone(D, bPrefix, [D.nchannels D.nsamples D.ntrials]); % save meeg with badchannel indices in Badchannel directory
Dbad(:,:,:) = D(:,:,:);
save(Dbad);

