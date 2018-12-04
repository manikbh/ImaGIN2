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

try
    bIdx = load(strcat(bPrefix,'_bIdx.txt'));
catch
    disp('Bad channel indices was not loaded');
    bIdx = [];
end
elec= sensors(D,'eeg');  % add channels without positions into bad channels
pos = elec.elecpos; 
chanLbs = elec.label;
idxNaN = find(isnan(pos(:,1)));
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
        fprintf(badchaFile, '%d %s\n', bIdx(i), chanLbs{bIdx(i)});
    end    
end
fclose(badchaFile);
tmpIdx = D.badchannels; % old badchannels indices
D = badchannels(D,tmpIdx,0); % reset meeg object bad indices

if ~isempty(bIdx)
    D = badchannels(D, bIdx,1); % add new badchannel indeces in meeg object
end
Dbad = clone(D, bPrefix, [D.nchannels D.nsamples D.ntrials]); % save meeg with badchannel indices in Badchannel directory
Dbad(:,:,:) = D(:,:,:);
save(Dbad);

