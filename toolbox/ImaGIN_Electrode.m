function ImaGIN_Electrode(S)
% Set electrode positions.

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
% Authors: Olivier David, Francois Tadel

% % Test findChannel()
% ListCSV = {'A01', 'A02', 'A03', 'A04', 'A05', 'A18', 'Bp01', 'Bp02', 'BP01', 'BP04', 'Pp01', 'T101', 'T111', 'V101', 'V201'};
% disp(['pp1: ',   ListCSV{findChannel('pp1', ListCSV)}])
% disp(['A01: ',   ListCSV{findChannel('A01', ListCSV)}])
% disp(['a01: ',   ListCSV{findChannel('a01', ListCSV)}])
% disp(['a1: ',    ListCSV{findChannel('a1', ListCSV)}])
% disp(['a 1: ',   ListCSV{findChannel('a 1', ListCSV)}])
% disp(['A18: ',   ListCSV{findChannel('A18', ListCSV)}])
% disp(['B''1: ',  ListCSV{findChannel('B''1', ListCSV)}])
% disp(['b''01: ', ListCSV{findChannel('b''01', ListCSV)}])
% disp(['Bp01: ',  ListCSV{findChannel('Bp01', ListCSV)}])
% disp(['Bp1: ',   ListCSV{findChannel('Bp1', ListCSV)}])
% disp(['bp1: ',   ListCSV{findChannel('bp1', ListCSV)}])
% disp(['B,1: ',   ListCSV{findChannel('B,1', ListCSV)}])
% disp(['B, 01: ',  ListCSV{findChannel('B, 01', ListCSV)}])
% disp(['BP01: ',  ListCSV{findChannel('BP01', ListCSV)}])
% disp(['Bp2: ',   ListCSV{findChannel('Bp2', ListCSV)}])
% disp(['B''2: ',   ListCSV{findChannel('B''2', ListCSV)}])
% disp(['B, 2: ',  ListCSV{findChannel('B, 2', ListCSV)}])
% disp(['BP4: ',   ListCSV{findChannel('BP4', ListCSV)}])
% disp(['T111: ',  ListCSV{findChannel('T111', ListCSV)}])
% disp(['T11: ',   ListCSV{findChannel('T11', ListCSV)}])
% disp(['V101: ',  ListCSV{findChannel('V101', ListCSV)}])
% disp(['V11: ',   ListCSV{findChannel('V11', ListCSV)}])
% disp(['V101: ',  ListCSV{findChannel('V101', ListCSV)}])
% disp(['V21: ',   ListCSV{findChannel('V21', ListCSV)}])
% return;

% Get file to edit
try
    t = S.Fname;
catch
    t = spm_select(Inf, '\.mat$', 'Select data file');
end
if isempty(t)
    return;
end
P = spm_str_manip(deblank(t(1,:)),'h');

% Read sensor names
Position = [];
try
    Name = S.Name;
catch
    try
        filename = S.filenameName;
    catch
        filename = spm_select(1, '\.(csv|txt)$', 'Select txt file for electrode names', {}, P);
    end
    % Test for input file
    if isempty(filename) || ~exist(filename, 'file')
        disp('ImaGIN> ERROR: Invalid input file name.');
        return
    end
    % Detect file type: _Name.txt/_Pos.txt or .csv
    [fPath, fBase, fExt] = fileparts(filename);
    % Read file
    switch lower(fExt)
        case '.txt'
            Name = readName(filename);
        case '.csv'
            [Name, Position] = readCsv(filename);
            if all(isnan(Position))
                warning('T1pre coordinates not found: possibly patient implanted before 2009.');
            end
        otherwise 
            disp('ImaGIN> ERROR: Invalid input file type.');
            return
    end
    if isempty(Name)
        disp('ImaGIN> ERROR: No contact names in input files.');
        return
    end
end

% Read sensor positions
if isempty(Position)
    try
        Position = S.Position;
    catch
        try 
            filename = S.filenamePos;
            Position = load(filename);
        catch
            filename = spm_select(1, '\.txt$', 'Select txt file for electrode positions', {}, P);
            Position = load(filename);
        end
    end
end


% Set positions
chNotFound = {};
chMatchLog = {};
for i0 = 1:size(t,1)
    T = deblank(t(i0,:));
    % Clone file if requested in input
    if (nargin >= 1) && isfield(S, 'FileOut') && ~isempty(S.FileOut)
        D = spm_eeg_load(T);
        D2 = clone(D, S.FileOut, [D.nchannels D.nsamples D.ntrials]); 
        D2(:,:,:) = D(:,:,:);
        save(D2);
        SpmFile = S.FileOut;
    else
        SpmFile = T;
    end
    
    % Load channel names
    SpmMat = load(SpmFile);
    Sensors = SpmMat.D.sensors.eeg;
    if length(unique(Sensors.label))~=length(Sensors.label)
        warning('Repeated label in the file.')
    end    
    % Loop on all channels available in the file
    for i1 = 1:length(Sensors.label)
        sensLtmp = Sensors.label{i1};
        sensLtmp(ismember(double(sensLtmp),[',' ';' '-' '_'])) ='';
        Sensors.label{i1} = sensLtmp;
        iChanPos = findChannel(Sensors.label{i1}, Name, 'all_upper');
        % If the channel was already found in the list before: check the best option based on the case
        if ~isempty(iChanPos) && ~isempty(chMatchLog)
            iPrevious = find(strcmp(Name{iChanPos}, chMatchLog(:,2)));
            if ~isempty(iPrevious) 
                % If the new channel has strictly the same case, or if it corresponds to a replaced "prime": remove the previous match
                if isequal(Sensors.label{i1}, Name{iChanPos}) || ...
                   (any(Sensors.label{i1} == '''') && strcmp(strrep(upper(Sensors.label{i1}), '''', 'p'), Name{iChanPos}))
                    disp(['ImaGIN> WARNING: Channel name conflict: ' Sensors.label{i1} ' matched with ' Name{iChanPos} ', ' chMatchLog{iPrevious,1} ' discarded.']);
                    chMatchLog(iPrevious,:) = [];                
                % If both labels are the same regardless of the case, raises an error.
                elseif strcmpi(chMatchLog{iPrevious,1},Sensors.label{i1})
                    error(['Duplicated label found: ' Sensors.label{i1}]);                                        
                else
                    disp(['ImaGIN> WARNING: Channel name conflict: ' chMatchLog{iPrevious,1} ' matched with ' Name{iChanPos} ', ' Sensors.label{i1} ' discarded.']);
                    iChanPos = [];
                end
            end
        end
        % If channel was found
        if ~isempty(iChanPos)
            % Copy position
            Sensors.elecpos(i1,:) = Position(iChanPos,:);
            Sensors.chanpos(i1,:) = Position(iChanPos,:);
            % Copy channel name from input name file (ADDED BY FT 5-Oct-2018)
            chMatchLog{end+1,1} = Sensors.label{i1};
            chMatchLog{end,2} = Name{iChanPos};
            Sensors.label{i1} = Name{iChanPos};
        else
            disp(['ImaGIN> WARNING: ' Sensors.label{i1} ' not assigned']);            
            Sensors.elecpos(i1,:) = NaN;
            Sensors.chanpos(i1,:) = NaN;
            chNotFound{end+1} = Sensors.label{i1};
        end
    end
    if isempty(chMatchLog)
        chMatchLog = repmat({''},1,2);        
    end
    % Electrodes present in the CSV but not found in the SEEG recordings
    [~, ~, chTagsCSV] = ImaGIN_select_channels(Name);
    [~, ~, chTagsSEEG] = ImaGIN_select_channels(chMatchLog(:,2)');
    elecUnused = setdiff(unique(chTagsCSV), unique(chTagsSEEG));
    
    % Add entries in a NEW log file
    if ~isempty(chMatchLog) 
        ImaGIN_save_log(SpmFile, 'Positions added for channels:', chMatchLog(:,1));
    end
    if ~isempty(chNotFound) 
        ImaGIN_save_log(SpmFile, 'Unmatched SEEG channels:', chNotFound);
        % save unmatched channels in separate .txt file
        [unPath, unFile, ~] = fileparts(SpmFile);
        if numel(chNotFound) > numel(cell2mat(regexpi(chNotFound,'ecg')))
            try
                unmatchedFileName  = fullfile(unPath, ['contactsunmatched_', unFile,'.txt']);
                fid = fopen(unmatchedFileName,'w');
                for i = 1:length(chNotFound)
                    if isempty(regexpi(chNotFound{i},'ecg'))
                        fprintf(fid,'%s\n', chNotFound{i});
                    end
                end
                fclose(fid);
            catch
                disp('Unmatched channels names SEEG-CSV not saved.')
            end
        end
    end
    if ~isempty(elecUnused) 
        ImaGIN_save_log(SpmFile, 'Unmatched CSV electrodes:', elecUnused);
    end
    
    % Match file: Correspondance between SEEG-CSV-LENA conventions
    try
        if isfield(S, 'FileTxtOut') && ~isempty(S.FileTxtOut)
            try
                fid = fopen(S.FileTxtOut,'w');
                fprintf(fid,'SEEG,CSV\n');
                for i = 1:size(chMatchLog,1)
                    fprintf(fid,'%s,%s\n', chMatchLog{i,1}, chMatchLog{i,2});
                end
                fclose(fid);
            catch
                disp('Log with matched channels names SEEG-CSV not saved.')
            end
        end
    end
    
    % Replace channel definition in input .mat/.dat file
    SpmMat.D.sensors.eeg = Sensors;
    
    % Replace labels in D.channels
    for iChan = 1:length(SpmMat.D.channels)
        spmLtmp = SpmMat.D.channels(iChan).label;
        spmLtmp(ismember(double(spmLtmp),[',' ';' '-' '_'])) ='';
        SpmMat.D.channels(iChan).label = spmLtmp;
        
        iChanMatch = find(strcmpi(SpmMat.D.channels(iChan).label, chMatchLog(:,1)));
        if (length(iChanMatch) == 1)
            SpmMat.D.channels(iChan).label = chMatchLog{iChanMatch,2};
        end
    end
    
    % Replace labels in events
    for iEvt = 1:length(SpmMat.D.trials.events)
        % Get channel name from event name
        EvtName = SpmMat.D.trials.events(iEvt).type;
        [chLabel1, chLabel2, noteNameNew, chInd1, chInd2] = ImaGIN_CleanEventName(EvtName);  
        if ~isempty(chInd1) && ~isempty(chInd2)
            if str2double(chInd1) <  str2double(chInd2)
                if str2double(chInd1) +1 ~= str2double(chInd2)
                    chInd2tmp = num2str(str2double(chInd1) + 1);
                    if numel(chInd1) ~= numel(chInd2)
                        chInd2tmp = strcat('0',chInd2tmp);
                    end
                    chLabel2 = strrep(chLabel2,chInd2,chInd2tmp);
                    chInd2   = chInd2tmp;
                end
            elseif str2double(chInd1) >  str2double(chInd2)
                if str2double(chInd1) ~= str2double(chInd2) + 1
                    chInd1tmp = num2str(str2double(chInd2) + 1);
                    if numel(chInd2) ~= numel(chInd1)
                        chInd1tmp = strcat('0',chInd1);
                    end
                    chLabel1 = strrep(chLabel1,chInd1,chInd1tmp);
                    chInd1   = chInd1tmp;
                end 
            end
        end
        % Replace with matching CSV name
        % This section now uses "Name" extracted from the csv to consider all electrode labels instead of the ones in "chMatchLog" as it used to do because it only considered electrodes which also recorded.
        % Note by A. Boyer - 13/02/2020
        elec_labels = Name';
        elec_labels_no_primes = strrep(elec_labels,'p',''''); 
        csv_all_electrodes = [elec_labels elec_labels_no_primes];      
        
        [iChanMatch1,~] = find(strcmpi(chLabel1, csv_all_electrodes));
        [iChanMatch2,~] = find(strcmpi(chLabel2, csv_all_electrodes));        
        if isempty(iChanMatch1) && ~isempty(chInd1)
            tmpchLabel1 = strrep(chLabel1, chInd1, num2str(str2double(chInd1)));
            [iChanMatch1,~] = find(strcmpi(tmpchLabel1, csv_all_electrodes));
            if isempty(iChanMatch1)
                tmpchLabel1 = strrep(chLabel1, chInd1, ['0' chInd1]);
                [iChanMatch1,~] = find(strcmpi(tmpchLabel1, csv_all_electrodes));
            end   
        end        
        if isempty(iChanMatch2) && ~isempty(chInd2)
            tmpchLabel2 = strrep(chLabel2, chInd2, num2str(str2double(chInd2)));
            [iChanMatch2,~] = find(strcmpi(tmpchLabel2, csv_all_electrodes));
            if isempty(iChanMatch2)
                tmpchLabel2 = strrep(chLabel2, chInd2, ['0' chInd2]);
                [iChanMatch2,~] = find(strcmpi(tmpchLabel2, csv_all_electrodes));
            end
        end
        
        if ~isempty(iChanMatch1)
            if sum(iChanMatch1==iChanMatch1(1)) == numel(iChanMatch1)
                iChanMatch1 = iChanMatch1(1);
            else
                warning('Same label found for multiple electrodes');
            end
        end
        if ~isempty(iChanMatch2)
            if sum(iChanMatch2==iChanMatch2(1)) == numel(iChanMatch2)
                iChanMatch2 = iChanMatch2(1);
            else
                warning('Same label found for multiple electrodes');
            end
        end

        if (length(iChanMatch1) == 1)
            noteNameNew = strrep(noteNameNew, chLabel1,csv_all_electrodes{iChanMatch1,1}); % chMatchLog{iChanMatch1,2}
        end
        if (length(iChanMatch2) == 1)
            noteNameNew = strrep(noteNameNew, chLabel2,csv_all_electrodes{iChanMatch2,1}); % chMatchLog{iChanMatch2,2}            
        end
        if (length(iChanMatch1) == 1) &&  (length(iChanMatch2) == 1)
        SpmMat.D.trials.events(iEvt).type = noteNameNew;
        end
    end
    if ~(exist('csv_all_electrodes','var') == 1)
        error('Could not match notes with channels found in the implantation file.')
    end  
    csv_struct.csv_labels = csv_all_electrodes(:,1);
    SpmMat.D.other = csv_struct; % Add an extra field to the .mat so we have a listing of all channel labels in the csv.
    % Update existing .mat file
    save(SpmFile, '-struct', 'SpmMat');
    save(spm_eeg_load(SpmFile)); % SPM's save to create a valid SPM object
end

end


%% Read name from file
function Name = readName(filename)
    fid = fopen(filename);
    i1 = 1;
    while 1
        tmp = fgetl(fid);
        if isequal(tmp, -1) || isempty(tmp)
            break;
        end
        Name{i1} = tmp;
        i1 = i1 + 1;
    end
    fclose(fid);
end


%% Read name and position from .csv
function [Name, Pos] = readCsv(filename)
    Name = {};
    Pos = [];
    % Open file
    fid = fopen(filename);
    i1 = 1;
    % Skip two header lines
    fgetl(fid);
    fgetl(fid);
    % Read column names (tab-separated file)
    colNames = strsplit(fgetl(fid), '\t');
    iColName = find(strcmpi(colNames, 'contact'));
    iColPos = find(strcmpi(colNames, 'T1pre Scanner Based'));
    if isempty(iColName) || isempty(iColPos)
        disp('ImaGIN> ERROR: Invalid CSV file: No columns "T1pre Scanner Based" or "contact".');
        return;
    end
    % Read contact by contact
    while 1
        % Read line
        tmp = fgetl(fid);
        if isequal(tmp, -1) || isempty(tmp)
            break;
        end
        % Split by columns (tab-separated file)
        tmp = strsplit(tmp, '\t');
        Name{i1} = tmp{iColName};
        if ~isempty(str2num(tmp{iColPos}))
            Pos(i1,1:3) = str2num(tmp{iColPos});
        else
            Pos(i1,1:3) = NaN;
        end
        i1 = i1 + 1;
    end
    fclose(fid);
end


%% Find channel name in a list
function iChanPos = findChannel(Label, List, caseType)
    % Test three different case versions
    % The goal is to be able to handle a List containing both electrodes L' ("Lp") and LP ("LP")
    % The List provided here comes from a .csv with this convention: only capitals except for the "p" for "'"
    if (nargin < 3) || isempty(caseType)
        iChanPos = findChannel(Label, List, 'no_change');
        if isempty(iChanPos)
            iChanPos = findChannel(Label, List, 'all_upper');
        end
        if isempty(iChanPos)
            iChanPos = findChannel(Label, List, 'upper_except_p');
        end
        return;
    end
    
    % Remove spaces
    Label(Label == ' ') = [];
    % Switch case type
    switch (caseType)
        case 'no_change'
            % Nothing to change
        case 'all_upper'
            Label = upper(Label);
        case 'upper_except_p'
            iP = find(Label == 'p');
            Label = upper(Label);
            if any(iP >= 2)
                Label(iP(end)) = 'p';
            end
    end
    % Replacing ' with p
    Label = strrep(Label, '''', 'p');
    List = strrep(List, '''', 'p');
    % Replacing , with p
    Label = strrep(Label, ',', 'p');
    List = strrep(List, ',', 'p');
    % Look for channel in position file
    iChanPos = find(strcmp(Label, List));
    if ~isempty(iChanPos)
        return;
    end    
    % Channel not found: try adding missing 0 (A1=>A01, T11=>T101)
    % Find the last letter in the name
    iLastLetter = find(~ismember(Label, '0123456789'), 1, 'last');
    % If there is not electrode label or no index: wrong naming
    if isempty(iLastLetter) || (iLastLetter == length(Label))
        return;
    end
    % Find channel index
    chLabel = Label(1:iLastLetter);
    chInd = str2num(Label(iLastLetter+1:end));
    % If label is > 100: Add 1 or 2 at the end of the label
    if (chInd > 220)
        return;
    elseif (chInd > 200)
        chLabel = [chLabel '2'];
        chInd = chInd - 200;
    elseif (chInd > 100)
        chLabel = [chLabel '1'];
        chInd = chInd - 100;
    elseif (chInd > 20)
        chLabel = [chLabel '2'];
        chInd = chInd - 20;
    end
    
    % Try with a zero
    iChanPos = find(strcmp(sprintf('%s%02d', chLabel, chInd), List));
    if ~isempty(iChanPos)
        return;
    end
    % Try without the zero
    iChanPos = find(strcmp(sprintf('%s%d', chLabel, chInd), List));
    if ~isempty(iChanPos)
        return;
    end

    % For channel indices between 11 and 19, try to interpret them as if electrode name was ending with a 1
    if (chInd >= 11) && (chInd <= 19)
        iChanPos = find(strcmp(sprintf('%s%02d', [chLabel '1'], chInd - 10), List));
        if ~isempty(iChanPos)
            return;
        end
    end
end
