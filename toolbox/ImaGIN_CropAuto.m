function D = ImaGIN_CropAuto(S)
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


% This function does not use separate csv and will delete stimulations done on
% electrodes (series of contacts) not recorded at all. Ideally, a csv with
% all existing contacts should be provided as input. Note by O. David on
% 19/12/19.
% The workaround introduced above modifies the naming convention of the
% files/stimulations which is very sensible and impacts the rest of the
% pipeline. Changes were reverted and modifications were made on previous
% steps so all the stimulations found in the notes are considered, even if
% one of the stimulating contact was never used for recording.
% Note by A. Boyer on 12/02/2020

if (nargin >= 1) && isfield(S, 'dataset') && ~isempty(S.dataset)
    sFile = S.dataset;
else
    sFile = spm_select(1, '\.mat$', 'Select EEG mat file');
    if isempty(sFile)
        return;
    end
end

if (nargin >= 1) && isfield(S, 'DirFileOut') && ~isempty(S.DirFileOut)
    DirOut= S.DirFileOut;
else
    [DirOut, ~, ~] = fileparts(sFile);
end

if (nargin >= 1) && isfield(S, 'StimName') && ~isempty(S.StimName)
thisN = S.StimName; % stim event to crop individually 
else
    thisN = '';
end

% method to be used in ImaGIN_StimDetect
try
    StimDetectVersion = S.StimDetectVersion ;
catch
    StimDetectVersion = 4 ;
end


D = spm_eeg_load(sFile); % Load the converted file .mat


%% Find existing crop files in DirOut
matExist = dir(fullfile(DirOut,'*.mat')); %Delete all cropped text file
matExist = {matExist.name};
matDelCount = 0;
matDeleted = {};
%%
evt = events(D);
evsize = size(evt,2);        % Number of events
Notes  = cell(1,evsize);     % Events labels
Time   = zeros(1,evsize);
totTime= max(time(D));
Time0  = min(time(D));
% Extract events properties (label and time in sampling)
for i = 1: evsize
    Notes{i}  = evt(i).type;
    Time(1,i) = evt(i).time;
end
nHz  = 5;  % Accepted stim frequency
minStim = 3; 
bgnTime = zeros(1,evsize); % Beginning of event
endTime = zeros(1,evsize); % End of event

%%
% Read notes to keep only those related to a stimulation
KeepEvent=[];
for c=1:evsize % Navigate all available events
    
    xpr1  = '\w*hz_\w*';
    xpr2  = '\w*stim\w*';
    xpr2b  = '\w*CONNECT TO:\w*';
    xpr3  = '\w*mA\w*';
    xpr4  = '\w*50.0hz\w*';
    xpr5  = '\w*50hz\w*';
    xpr6  = '\w*50 hz\w*';
    
    xpr4b  = '\w*55.0hz\w*';
    xpr5b  = '\w*55hz\w*';
    xpr6b  = '\w*55 hz\w*';
    % present in '/gin/data/database/02-raw/0007BUC/SEEG/StimLF_2013-11-13/Electrodes/electrodes_0007BUC_Raw_2.mat' ; 'P10-P11_2mA_1011Hz_3000us'
    % while in fact, it is 1Hz
    % and using 1011Hz makes a poor detection 
    xpr7a =  '\w*Yf\w*';
    xpr7  = '\w*alarme\w*';
    xpr8  = '\w*SE1Hz\w*';
    xpr9  = '\w*SE 1Hz\w*';
    xpr10  = 'crise';
    xpr11  = '\w*OFF\w*';
    xpr12  = '\w*t65535\w*';
    xpr12b = '\w*Test\w*';
    xpr13  = '\w*Stopw*';
    
    [ds,di]=regexp(Notes{c},'\d*','Match');
    if ~isempty(di)
        if ~isempty(regexpi(Notes{c},xpr4))       || ~isempty(regexpi(Notes{c},xpr5)) || ...
                ~isempty(regexpi(Notes{c},xpr6))  || ~isempty(regexpi(Notes{c},xpr7)) || ~isempty(regexpi(Notes{c},xpr7a)) ||...
                ~isempty(regexpi(Notes{c},xpr8))  || ~isempty(regexpi(Notes{c},xpr9)) || ...
                ~isempty(regexpi(Notes{c},xpr12)) || ~isempty(regexpi(Notes{c},xpr4b))|| ...
                ~isempty(regexpi(Notes{c},xpr12b))|| ~isempty(regexpi(Notes{c},xpr13))|| ...
                ~isempty(regexpi(Notes{c},xpr5b)) || ~isempty(regexpi(Notes{c},xpr6b))|| ...
                ~isempty(regexpi(Notes{c},xpr11)) || strcmpi(Notes{c}(1:min([length(Notes{c}) 5])),xpr10)
        elseif ~isempty(regexpi(Notes{c},xpr1))
            KeepEvent=[KeepEvent c];
        elseif ~isempty(regexpi(Notes{c},xpr2))
            KeepEvent=[KeepEvent c];
        elseif ~isempty(regexpi(Notes{c},xpr2b))
            KeepEvent=[KeepEvent c];
        elseif ~isempty(regexp(Notes{c},xpr3,'ONCE'))
            KeepEvent=[KeepEvent c];
        elseif ismember(lower(Notes{c}(1)),['a':'z']) && di(1)<=4 && ~strcmp(regexprep(Notes{c},' ',''),'SE1Hz') && ~strcmp(regexprep(Notes{c},' ',''),'SE50Hz')
            KeepEvent=[KeepEvent c];
        end
    end
end

% Clean KeepEvent (Remove stim over frequency max, and in case a note was written for every pulse, we keep only the first of the series)
IndexToRemove=[];
for i1=1:length(KeepEvent)-1
    if strcmp(Notes(KeepEvent(i1)),Notes(KeepEvent(i1+1)))
        IndexToRemove=[IndexToRemove i1+1];
    end
end
for i1=1:length(KeepEvent)
    
    noteName = strrep(char(Notes{KeepEvent(i1)}), ' ','_');
    noteName = strrep(noteName,'�sec','us');
    noteName = strrep(noteName,'�s','us');
    noteName(~ismember(double(noteName),['A':'Z' 'a':'z' '_' '.' '''' 'µ' '-' '0':'9'])) ='_';
    noteName = strrep(noteName,'_+','_'); noteName = strrep(noteName,'µ','u');
    noteName = strrep(noteName,'usec','us');
    noteName = strrep(noteName,'Stim_Start_',''); %YUQ & MIL notes
    noteName = strrep(noteName,'Stim_Stop_','');  %YUQ notes
    noteName = strrep(noteName,'-','_'); noteName = strrep(noteName,'__','_');
    noteName = strrep(noteName,',',''); noteName = strrep(noteName,'_mA_','_');
    noteName = strrep(noteName,'sec','us'); noteName = strrep(noteName,'_us','us');
    noteName = strrep(noteName,'AA','A'); noteName = strrep(noteName,'_MA_','_'); %some MIL notes
    noteName = strrep(noteName,'stim',''); noteName = strrep(noteName,'Stim','');
    noteName = strrep(noteName,'CONNECT TO:',''); noteName = strtrim(noteName); % for MNI notes
    noteName = strrep(noteName,'TextNote:',''); % for BRN datasets
    noteName = strrep(noteName,'Stop','');    noteName = strrep(noteName,'Start','');
    noteName = strrep(noteName,'1011Hz','1Hz'); noteName = strrep(noteName,'uus','us');
    noteName = strrep(noteName,'Impulsions_biphasiques',''); noteName = strrep(noteName,'CDuruue','');
    
    [ds,di] = regexp(noteName,'\d*','Match');
    xsub0 = noteName(1:(di(1)+numel(ds{1})-1));
    rxp1  = '[-+]?(\d*[.])?\d+mA'; rxp2  = '[-+]?(\d*[.])?\d+Hz';
    rxp3  = '[-+]?(\d*[.])?\d+us'; rxp4  = '[-+]?(\d*[.])?\d+s';
    xsub1 = regexpi(noteName,rxp1,'match');
    if isempty(xsub1), xsub1 = '0mA';
        xsub1 = cellstr(xsub1);
    end
    xsub2 = regexp(noteName,rxp2,'match');
    if isempty(xsub2) %|| strcmp(xsub2,'0Hz')
        xsub2 = '0Hz';
        xsub2 = cellstr(xsub2);
        S.StimFreq = 1;
        StimFreq = 1;
    elseif strcmp(xsub2,'0Hz')  %we assume 0Hz should read 0.26 (FRE)
        xsub2 = '026Hz';
        xsub2 = cellstr(xsub2);
        S.StimFreq = 0.26;
        StimFreq = 0.26;
    else
        % We assume a leading 0 implies a stimulation frequency <1 (e.g. 026 -> 0.26).
        % Only used for the detection of stimulations (ImaGIN_StimDetect.m) Boyer.A 25/08/2020
        freq_string = xsub2{1}(1:end-2);
        match_freq = regexp(freq_string,'^0+(?<value>[0-9]+)$','once','names');
        if ~isempty(match_freq)
            StimFreq=str2double(['0.' match_freq.value]);
        else
            StimFreq=str2double(freq_string);
        end
    end
    if StimFreq>nHz
        IndexToRemove=[IndexToRemove i1];
    end
end
KeepEvent=KeepEvent(setdiff(1:length(KeepEvent),IndexToRemove));



%% ------------------------------------------------
% Case only specific stim events are to be cropped
if ~isempty(thisN)
    idxN = strcmp(Notes,thisN);
    idxN= find(idxN);
    KeepEvent = find(KeepEvent==idxN);
end
%% ------------------------------------------------
for c = 1:length(KeepEvent) % Navigate all stim events
    
    % Detect stimulations and stimulation indices & save in text file
    [pth,matFile,~]= fileparts(sFile);
    clear S
    S.Fname = fullfile(pth, matFile);
    S.Channels = [];
    S.StimStart= Time(1,KeepEvent(c))-2;
    if c < length(KeepEvent)
        thisEnd = min([Time(1,KeepEvent(c))+180 Time(1,KeepEvent(c+1))-0.5]);
        if thisEnd <= totTime
            S.StimEnd = thisEnd;
        else
            continue;
        end
    else
        thisEnd = min([totTime Time(1,KeepEvent(c)) + 120]);
        if thisEnd <= totTime && S.StimStart < thisEnd
            S.StimEnd  = thisEnd;
        else 
            continue
        end
    end
    S.StimContinuous= 0;
    S.EvtName  = Notes{KeepEvent(c)};
    noteName = strrep(char(Notes{KeepEvent(c)}), ' ','_');
    noteName = strrep(noteName,'�sec','us');
    noteName = strrep(noteName,'�s','us');
    noteName(~ismember(double(noteName),['A':'Z' 'a':'z' '_' '.' '''' 'µ' '-' '0':'9'])) ='_';
    noteName = strrep(noteName,'_+','_'); noteName = strrep(noteName,'µ','u');
    noteName = strrep(noteName,'usec','us'); 
    noteName = strrep(noteName,'Stim_Start_',''); %YUQ & MIL notes
    noteName = strrep(noteName,'Stim_Stop_','');  %YUQ notes    
    noteName = strrep(noteName,'-','_'); noteName = strrep(noteName,'__','_');    
    noteName = strrep(noteName,',',''); noteName = strrep(noteName,'_mA_','_');
    noteName = strrep(noteName,'sec','us'); noteName = strrep(noteName,'_us','us');
    noteName = strrep(noteName,'AA','A'); noteName = strrep(noteName,'_MA_','_'); %some MIL notes 
    noteName = strrep(noteName,'stim',''); noteName = strrep(noteName,'Stim',''); 
    noteName = strrep(noteName,'CONNECT TO:',''); noteName = strtrim(noteName); % for MNI notes
    noteName = strrep(noteName,'TextNote:',''); % for BRN datasets 
    noteName = strrep(noteName,'Stop','');    noteName = strrep(noteName,'Start','');
    noteName = strrep(noteName,'1011Hz','1Hz'); noteName = strrep(noteName,'uus','us');
    noteName = strrep(noteName,'Impulsions_biphasiques',''); noteName = strrep(noteName,'CDuruue','');
    [numZ, numZI] = regexp(noteName,'\d*','Match');
    keepN = '';
    try
        fundc = strfind(noteName,'_');
        lNumb = strfind(noteName,noteName(1:fundc(1)-1));
        iLastLetter = find(~ismember(noteName(1:fundc(1)-1),   '0123456789'), 1, 'last');
        keepN = noteName(1:fundc(1)-1);
        if(numel(lNumb)) == 2 && ~strcmp(keepN,'A') && ~strcmp(keepN,'H') && ~isempty(keepN)            
            noteName = strrep(noteName,keepN,'CHNAME'); 
        elseif iLastLetter > 0 && iLastLetter > numZI(1)
             keepN = noteName(1:iLastLetter);
             noteName = strrep(noteName,keepN,'CHNAME');
        end  
    end   

    noteName = strrep(noteName,'.0',''); noteName = strrep(noteName,'.','');
    idScore = strfind(noteName,'_');
    if ~isempty(idScore)
        if idScore(1) < numZI(1)
           noteName(idScore(1)) = '';
        end
    end
    %% check if stim electr numbers are concatenated without space or -
    
    numbr = regexp(noteName,'\d*','Match');
    if numel(numbr) > 2 && numel(numbr) <= 5
        if str2double(numbr(2))~= str2double(numbr(1)) + 1 && str2double(numbr(1))~= str2double(numbr(2)) + 1 
            if numel(numbr{1}) == 2
                 if numel(numbr{2}) ~= 3 && numel(numbr{2}) ~= 2 && numel(numbr{2}) ~= 4
                     if str2double(numbr{1}(1)) + 1 == str2double(numbr{1}(2)) || str2double(numbr{1}(1)) == str2double(numbr{1}(2)) + 1
                         elecno = strcat(numbr{1}(1),'_',numbr{1}(2));
                         noteName = strrep(noteName,numbr{1},elecno);
                     end
                elseif numel(numbr{2}) == 3                      
                    if str2double(numbr{2}(1:2)) == str2double(numbr(1)) + 1 || str2double(numbr{2}(1:2)) + 1 == str2double(numbr(1)) 
                        elecno = strcat(numbr{2}(1:2),'_',numbr{2}(3));
                        noteName = strrep(noteName,numbr{2},elecno);
                    elseif str2double(numbr{2}(1)) == str2double(numbr(1)) + 1 || str2double(numbr{2}(1)) + 1 == str2double(numbr(1)) 
                        elecno = strcat(numbr{2}(1),'_',numbr{2}(2:3));
                        noteName = strrep(noteName,numbr{2},elecno);                        
                    end
                 elseif numel(numbr{2}) == 4
                     if str2double(numbr{2}(1:2)) == str2double(numbr(1)) + 1 || str2double(numbr{2}(1:2)) + 1 == str2double(numbr(1))
                         elecno = strcat(numbr{2}(1:2),'_',numbr{2}(3:4));
                         noteName = strrep(noteName,numbr{2},elecno);
                     end
                 elseif numel(numbr{2})== 2
                   if str2double(numbr{2}(1)) == str2double(numbr(1)) + 1 || str2double(numbr{2}(1)) + 1 == str2double(numbr(1))
                      elecno = strcat(numbr{2}(1),'_',numbr{2}(2)); 
                      noteName = strrep(noteName,numbr{2},elecno);
                   end
                end
            elseif numel(numbr{1}) == 3
                if str2double(numbr{1}(1))== 9 && str2double(numbr{1}(2:3))== 10 
                    elecno = strcat(numbr{1}(1),'_',numbr{1}(2:3));
                    noteName = strrep(noteName,numbr{1},elecno);
                elseif str2double(numbr{1}(1:2))== 10 && str2double(numbr{1}(2:3))== 9
                    elecno = strcat(numbr{1}(1:2),'_',numbr{1}(3));
                    noteName = strrep(noteName,numbr{1},elecno);
                elseif str2double(numbr{1}(1)) + 1 == str2double(numbr{1}(2)) || str2double(numbr{1}(1)) == str2double(numbr{1}(2)) + 1
                    elecno = strcat(numbr{1}(1),'_',numbr{1}(2),'_',numbr{1}(3));
                    noteName = strrep(noteName,numbr{1},elecno);
                end
            elseif numel(numbr{1}) == 4
                if str2double(numbr{1}(1:2)) + 1 ==  str2double(numbr{1}(3:4)) || str2double(numbr{1}(1:2)) ==  str2double(numbr{1}(3:4)) + 1
                    elecno = strcat(numbr{1}(1:2),'_',numbr{1}(3:4));
                    noteName = strrep(noteName,numbr{1},elecno);
                elseif str2double(numbr{1}(1))== 9 && str2double(numbr{1}(2:3)) == 10 
                    elecno = strcat(numbr{1}(1),'_',numbr{1}(2:3),'_',numbr{1}(4));
                    noteName = strrep(noteName,numbr{1},elecno);
                elseif str2double(numbr{1}(1:2))== 10 && str2double(numbr{1}(3)) == 9
                    elecno = strcat(numbr{1}(1:2),'_',numbr{1}(3),'_',numbr{1}(4));
                    noteName = strrep(noteName,numbr{1},elecno);
                end
            elseif numel(numbr{1}) == 5
                 if str2double(numbr{1}(1:2)) + 1 ==  str2double(numbr{1}(3:4)) || str2double(numbr{1}(1:2)) ==  str2double(numbr{1}(3:4)) + 1
                    elecno = strcat(numbr{1}(1:2),'_',numbr{1}(3:4),'_',numbr{1}(5));
                    noteName = strrep(noteName,numbr{1},elecno);
                 end
            elseif numel(numbr{1}) == 1 
                if numel(numbr{2}) == 2
                    if str2double(numbr{2}(1)) == str2double(numbr(1)) + 1
                        elecno = strcat(numbr{2}(1),'_',numbr{2}(2));
                        noteName = strrep(noteName,numbr{2},elecno);
                    else
                       elecno = strcat(numbr{2}(1),'_',numbr{2}(2));
                       noteName = strrep(noteName,numbr{2},elecno);
                    end
                elseif numel(numbr{2}) == 3                    
                    if str2double(numbr{2}(1)) == str2double(numbr(1)) + 1 || str2double(numbr{2}(1)) + 1 == str2double(numbr(1))
                        elecno = strcat(numbr{2}(1),'_',numbr{2}(2:3));
                        noteName = strrep(noteName,numbr{2},elecno);
                    else                        
                        if str2double(numbr(1)) + 1 == str2double(numbr{2}(1:2))
                            elecno = strcat(numbr{2}(1:2),'_',numbr{2}(3));
                            noteName = strrep(noteName,numbr{2},elecno);
                        end
                    end
                end
            end
        end
    elseif numel(numbr) == 2
        if ~isempty(strfind(noteName, [numbr{2} 'us'])) || ~isempty(strfind(noteName, [numbr{2} 'mA'])) ...
                || ~isempty(strfind(noteName, [numbr{2} 'Hz']))
            if numel(numbr{1}) == 2
                elecno = strcat(numbr{1}(1),'_',numbr{1}(2));
                if str2double(numbr{1}(1)) + 1 == str2double(numbr{1}(2)) || str2double(numbr{1}(1)) == str2double(numbr{1}(2)) + 1
                    noteName = strrep(noteName,numbr{1},elecno);
                end
            elseif ~isempty(strfind(numbr{1},'10'))  
                 noteName = strrep(noteName,'10','_10_');
            elseif numel(numbr{1}) == 4
                elecno = strcat(numbr{1}(1:2),'_',numbr{1}(3));
                noteName = strrep(noteName,numbr{1},elecno);
            elseif ~isempty(strfind(noteName, [numbr{1} 'mA'])) || ~isempty(strfind(noteName, [numbr{1} 'Hz'])) 
                if isempty(strfind(numbr{1},'.'))
                    if numel(numbr{1}) == 3
                        elecno = strcat(numbr{1}(1),'_',numbr{1}(2),'_',numbr{1}(3));
                        noteName = strrep(noteName,numbr{1},elecno);
                    elseif numel(numbr{1}) == 5
                        elecno = strcat(numbr{1}(1:2),'_',numbr{1}(3:4),'_',numbr{1}(5));
                        noteName = strrep(noteName,numbr{1},elecno);
                    end
                end
            end
        end
    end


  %%   
    numb = regexp(noteName,'\d*','Match');
    
    if numel(numb) >= 2
        idx1 = strfind(noteName,numb{1});
        idx2 = strfind(noteName,numb{2});
        cnbre1= numel(numb{1});
        cnbre2= numel(numb{2});
        if str2double(numb(1)) == str2double(numb(2)) - 1
            sfix = noteName(1:(idx1(1)+cnbre1)-1);
            noteName = strcat((sfix),noteName(idx2(1):end));
        elseif str2double(numb(1)) == str2double(numb(2)) + 1
            sfix = noteName(1:idx1(1)-1);
            noteName = strcat((sfix),num2str(numb{2}), num2str(numb{1}),noteName(idx2(1)+cnbre2:end));
        else
            if str2double(numb(1)) > str2double(numb(2))
                sfix = noteName(1:idx1(1)-1);
                noteName = strcat((sfix),num2str(numb{2}), num2str(numb{1}),noteName(idx2(1)+cnbre2:end));
            elseif str2double(numb(1)) < str2double(numb(2))
                sfix = noteName(1:idx1(1)-1);
                noteName = strcat((sfix),num2str(numb{1}), num2str(numb{2}),noteName(idx2(1)+cnbre2:end));
            end
        end

        xpr1  = '\w*Hz_\w*'; xpr2 = '\w*us_\w*';
        xpri1 = regexpi(noteName,xpr1); xpri2 = regexpi(noteName,xpr2);
    end
    noteName = strrep(noteName,'-','');  noteName = strrep(noteName,'_mA','mA');
    
    %% build .mat/.dat name
    
    try
        [ds,di] = regexp(noteName,'\d*','Match');
        xsub0 = noteName(1:(di(1)+numel(ds{1})-1));
        rxp1  = '[-+]?(\d*[.])?\d+mA'; rxp2  = '[-+]?(\d*[.])?\d+Hz';
        rxp3  = '[-+]?(\d*[.])?\d+us'; rxp4  = '[-+]?(\d*[.])?\d+s';
        xsub1 = regexpi(noteName,rxp1,'match'); 
        if isempty(xsub1), xsub1 = '0mA'; 
            xsub1 = cellstr(xsub1);
        end
        xsub2 = regexp(noteName,rxp2,'match'); 
        if isempty(xsub2) %|| strcmp(xsub2,'0Hz')
            xsub2 = '0Hz';
            xsub2 = cellstr(xsub2);
            S.StimFreq = 1;
            StimFreq = 1;
        elseif strcmp(xsub2,'0Hz')  %we assume 0Hz should read 0.26 (FRE)
            xsub2 = '026Hz';
            xsub2 = cellstr(xsub2);
            S.StimFreq = 0.26;
            StimFreq = 0.26;
        else           
            % We assume a leading 0 implies a stimulation frequency <1 (e.g. 026 -> 0.26).
            % Only used for the detection of stimulations (ImaGIN_StimDetect.m) Boyer.A 25/08/2020
            freq_string = xsub2{1}(1:end-2);
            match_freq = regexp(freq_string,'^0+(?<value>[0-9]+)$','once','names');
            if ~isempty(match_freq)
                StimFreq=str2double(['0.' match_freq.value]);
            else
                StimFreq=str2double(freq_string);
            end            
            S.StimFreq = StimFreq;
        end
        xsub3 = regexp(noteName,rxp3,'match');
        xsub4 = regexp(noteName,rxp4,'match');
        if ~isempty(xsub3)&& isempty(xsub4)
            xsub3 = cellstr(xsub3);
        elseif isempty(xsub3)&& ~isempty(xsub4)
            xsub4 = char(xsub4);
            if numel(xsub4) <= 3
                buff = num2str(str2double(xsub4)*1000);
                if numel(buff) <= 4
                    xsub4 = buff;
                end
            end
            xsub3 = cellstr(strcat(xsub4(1:end-1),'us'));
        else
            xsub3 = '0us';
        end
        
        FullN = strcat((xsub0),'_',xsub1,'_',xsub2,'_',xsub3);
        FullN = char(unique(FullN));
        FullN = FullN(1,:);
        if isempty(FullN)
            FullN = noteName;
        end
    catch
        FullN = noteName;
        S.StimFreq = 1;
        StimFreq = S.StimFreq;
    end
    
    noteName = strrep(FullN,'.',',');
    
    %%
    numb2 = regexp(noteName,'\d*','Match');
    if numel(numb2) >= 1
        idxn1 = strfind(noteName,numb2{1});          %OD
        subn1 = strrep(noteName(1:idxn1),'_','');
        noteName = char(strcat(subn1,noteName(idxn1+1:end)));
    end 
    
    noteName = strrep(noteName,'CHNAME',(keepN));  
    
    ptrn = ',';
    if strncmp(noteName,ptrn,1)
        noteName = char(noteName(2:end));
    end
    %%
    [~,tmpdi]=regexp(noteName,'\d*','Match');
    noteNameNew=noteName;
    % noteNameNew(1:tmpdi(1)-1)=upper(noteNameNew(1:tmpdi(1)-1)); 
    noteNameNew = strrep(noteNameNew,'''','p');     %to avoid ' in the name
        
    % Change of naming convention: Ap12_3mA_1Hz_1000us_1 becomes Ap1-Ap2_3mA_1Hz_1000us_1
    idxScore = strfind(noteNameNew,'_');
    if ~isempty(idxScore) && idxScore(1) > 1
        Label       = noteNameNew(1:idxScore(1)-1);
        iLastLetter = find(~ismember(Label, '0123456789'), 1, 'last');
    else
        Label       = '';
        iLastLetter = '';
    end
    if ~isempty(iLastLetter) && length(Label) >= iLastLetter
        chLabel =  Label(1:iLastLetter);
        chInd   =  Label(iLastLetter+1:end);
    else
        chLabel = '';
        chInd   = '';
    end
    
    if numel(chInd) == 2   
        chInd1 = chInd(1);
        chInd2 = chInd(2);
        chLabel1 = strcat(chLabel, chInd1);
        chLabel2 = strcat(chLabel, chInd2);        
        noteNameNew = strcat(chLabel1, '-',chLabel2, noteNameNew(idxScore(1):end));
    elseif numel(chInd) == 3        
        chInd1 = chInd(1);
        chLabel1 = strcat(chLabel, chInd1);        
        chInd2 = chInd(2:3);
        chLabel2 = strcat(chLabel,  chInd2);
        noteNameNew = strcat(chLabel1,'-',chLabel2, noteNameNew(idxScore(1):end));
    elseif numel(chInd) == 4
        chInd1 = chInd(1:2);
        chInd2 = chInd(3:4);
        chLabel1 = strcat(chLabel,chInd1);
        chLabel2 = strcat(chLabel, chInd2);
        noteNameNew = strcat(chLabel1, '-',chLabel2, noteNameNew(idxScore(1):end));
    elseif numel(chInd) == 6
        chInd1 = chInd(1:3);
        chInd2 = chInd(4:6);
        chLabel1 = strcat(chLabel, chInd1);
        chLabel2 = strcat(chLabel, chInd2);
        noteNameNew = strcat(chLabel1, '-',chLabel2, noteNameNew(idxScore(1):end));
    else
        chInd1 = '';
        chInd2 = '';
        chLabel1 = '';
        chLabel2 = '';
    end
    
    S.FileOut=  fullfile(DirOut, strcat(noteNameNew,'.txt'));   

    mat_file = load(sFile);
    csv_chanlabels = mat_file.D.other.csv_labels; % This new field is generated during the Electrode step so the .mat stores the channel labels found in the csv 
    
    iChanMatch1 = find(strcmpi(chLabel1, csv_chanlabels));
    iChanMatch2 = find(strcmpi(chLabel2, csv_chanlabels));
    
    if isempty(iChanMatch1) && ~isempty(chInd1)
        tmpchLabel1 = strrep(chLabel1, chInd1, num2str(str2double(chInd1)));
        iChanMatch1 = find(strcmpi(tmpchLabel1, csv_chanlabels));
    end
    if isempty(iChanMatch1) && ~isempty(chInd1)
        tmpchLabel1 = strrep(chLabel1, chInd1, ['_' chInd1]);
        iChanMatch1 = find(strcmpi(tmpchLabel1, csv_chanlabels));
    end  
    if isempty(iChanMatch2) && ~isempty(chInd2)
        tmpchLabel2 = strrep(chLabel2, chInd2, num2str(str2double(chInd2)));
        iChanMatch2 = find(strcmpi(tmpchLabel2, csv_chanlabels));
    end
    if isempty(iChanMatch2) && ~isempty(chInd2)
        tmpchLabel2 = strrep(chLabel2, chInd2, ['_' chInd2]);
        iChanMatch2 = find(strcmpi(tmpchLabel2, csv_chanlabels));
    end   
    
    %{  
    %% uncomment this section for some datasets
    % stimN1 = str2double(regexp(sfix,'\d*','Match'));
    % xstrN2 = strrep(xsub0,num2str(stimN1),'');
    % stimN2 = str2double(regexp(xstrN2,'\d*','Match'));
    % sfix2  = strrep(sfix,num2str(stimN1), num2str(stimN2));
    % selCh1 = find(strcmpi(elec.label,sfix));
    % selCh2 = find(strcmpi(elec.label,sfix2));
    % S.Channels = [selCh1,selCh2]; 
    % S.FindBadChannels=0;
    %}
    
    
    disp(KeepEvent(c)), disp(S.EvtName);
    S.StimDetectVersion = StimDetectVersion ;
    S.minStim= minStim;
    [stimTime,~,~] = ImaGIN_StimDetect(S);
    
    stimFq = StimFreq;
    
    % if numel(stimTime) >= seuilHz         %OD
    if numel(stimTime) >= minStim
        
        [~, stimTimeOut, ~] = fileparts(S.FileOut);
        
        %         if StimulationFreqU>20  %OD ugly fix to try to  detect 50 Hz stimulation
        %             stimFq = round(StimulationFreqU);
        %         else
        stimStep = stimTime(2:end) - stimTime(1:end-1);
        stimFq   = 1/median(stimStep); % median frequency within event : %OD corrected to median
        if stimFq>0.95
            stimFq=round(stimFq);
        end
        %         end
        
        % --------------------------------------------------
        % Define a file name for a dissected seeg
        rxp2  = '[-+]?(\d*[.])?\d+Hz';
        FileName = strcat(stimTimeOut,'_1.mat');
        stimHz = regexp(FileName,rxp2,'match');
        strFq = strcat(num2str(stimFq),'Hz');
        if stimFq>=1 && ~isempty(stimHz)
            if ~strcmp(stimHz{1},strFq) && ~strcmp(strFq,'0Hz')
                FileName = char(strrep(FileName,stimHz,strFq));
            end
        end
        % Check if the file exists i.e repeated stim
        if exist(fullfile(DirOut, FileName), 'file') ~= 2
            myFileName = FileName;
        else
            numstr = dir(fullfile(DirOut,[FileName(1:end-5),'*','.mat']));
            num    = size(numstr,1) + 1;
            myFileName = sprintf('%s%d.mat',FileName(1:end-5),num);
        end
        namePrefix = myFileName(1:end-4);
        % ----------------------------------------
        % Call ImaGIN_Crop
        bgnTime(KeepEvent(c)) = stimTime(1);
        endTime(KeepEvent(c)) = stimTime(end);
        
        % Select only the stimulations
        if stimFq <= nHz && numel(stimTime) >= minStim && numel(numb) >= 1 ...
                && ~strcmp(xsub2(1),'50Hz') && endTime(KeepEvent(c)) > Time(KeepEvent(c))
            clear S
            S.Job   = 'Manual';
            S.Fname = fullfile(pth, matFile);
            S.EventStart= Notes{KeepEvent(c)}; % Note n
            S.numbStart = KeepEvent(c); %%%
            S.numbEnd = KeepEvent(c); %%%
            if c > 1
                if endTime(KeepEvent(c-1)) == 0
                    endTime(KeepEvent(c-1)) = Time(KeepEvent(c-1));
                end
                if Time(KeepEvent(c)) - endTime(KeepEvent(c-1))-1 >= 40
                    setStart = 40;
                else
                    setStart = Time(KeepEvent(c)) - endTime(KeepEvent(c-1))-1;
                    if setStart<1
                        setStart=1;
                    end
                end
            elseif c == 1
                if Time(KeepEvent(c))-Time0 >= 40
                    setStart = 40;
                else
                    setStart = Time(KeepEvent(c)) - Time0;
                end
            end
            S.OffsetStart= setStart;
            S.EventEnd   = Notes{KeepEvent(c)}; % Note n
            S.OffsetEnd  = endTime(KeepEvent(c)) - Time(KeepEvent(c)) + 3;
            S.FileOut     = fullfile(DirOut, namePrefix);
            ImaGIN_Crop(S);
        end
    end
    
    if stimFq <= nHz && numel(stimTime) >= minStim && numel(numb) >= 1 ...
            && ~strcmp(xsub2(1),'50Hz') && endTime(KeepEvent(c)) > Time(KeepEvent(c))
        
        % Add events
        clear S2
        S2.Fname   = S.FileOut;
        S2.FileOut = S2.Fname;
        S2.Action  = 'Add';
        S2.Nevent  = 1;
        S2.EventName{1} = 'Stim';
        S2.EventFileName{1} = fullfile(DirOut, strcat(stimTimeOut,'.txt')); % txt file, time in seconds;
        ImaGIN_Events(S2);
        
        % Time origin
        clear S2
        S2.Fname   = S.FileOut;
        S2.FileOut = S2.Fname;
        S2.EventRef= 'Stim';
        S2.Offset  = 0;
        D = ImaGIN_TimeZero(S2);
        cropFileNameMat = fullfile(DirOut, strcat(namePrefix,'.mat'));
        cropFileNameDat = fullfile(DirOut, strcat(namePrefix,'.dat'));
        scoreIDX = strfind(namePrefix,'_');
        if numel(scoreIDX) == 4
            cropNumber = str2double(namePrefix(scoreIDX(4)+1:end));
            if ~isnan(cropNumber) && cropNumber > 1
                
                if exist(cropFileNameMat,'file')== 2
                    Dcrop = spm_eeg_load(cropFileNameMat);
                    refCropFileNameMats = dir(fullfile(DirOut, strcat(namePrefix(1:scoreIDX(4)),'*.mat')));
                    for i = 1:size(refCropFileNameMats,1)
                        refCropFileNameMat  = fullfile(DirOut, refCropFileNameMats(i).name);
                        if ~strcmp(refCropFileNameMat,cropFileNameMat)
                            if exist(refCropFileNameMat,'file')== 2 && exist(cropFileNameMat,'file')== 2
                                Dref  = spm_eeg_load(refCropFileNameMat);
                                if isequal(Dref(:,:,:),Dcrop(:,:,:))
                                    delete(cropFileNameMat);
                                    delete(cropFileNameDat);
                                end
                                clear Dref;
                            end
                        end
                    end
                end
            end
        end
        
        if exist(cropFileNameMat, 'file') == 2 && (isempty(iChanMatch1) || isempty(iChanMatch2))
            matDelCount =  matDelCount + 1;
            matDeleted(end+1) = {namePrefix};
        elseif exist(cropFileNameMat, 'file') == 2 && (str2double(chInd2) ~= str2double(chInd1) + 1)
            matDelCount =  matDelCount + 1;
            matDeleted(end+1) = {namePrefix};
        end
    end
end
%%
clc;
f = dir(fullfile(DirOut,'*.mat')); %Delete all cropped text file and .mat file with no existing stim channels
f = {f.name};
for k=1:numel(f)
    [~, matFileName, ~] = fileparts(f{k});
    nbrSc = strfind(matFileName,'_');
    if ~isempty(nbrSc)
        txtFileName = matFileName(1:nbrSc(end)-1);
        if exist(fullfile(DirOut,strcat(txtFileName,'.txt')),'file')== 2
            delete(fullfile(DirOut,strcat(txtFileName,'.txt')));
        end
    end
end
nf = dir(fullfile(DirOut,'*0us.txt'));
nf = {nf.name};

realStims = length(KeepEvent);
realCrps  = numel(f)-numel(matExist);
stimTxt   = numel(nf) + realCrps + matDelCount;
misCrop = 0;
nanStim = 0;

if ~isempty(nf)
    for k2 = 1:numel(nf)
        fID = fopen(fullfile(DirOut, nf{k2}));
        res = {};
        while ~feof(fID)
            res{end+1,1} =fgetl(fID);
        end
        fclose(fID);
        nbrLines = numel(res);
        if nbrLines > minStim
            misCrop = misCrop + 1;
            fprintf('%s not cropped\n',nf{k2});
        else
            nanStim = nanStim +1;
            fprintf('%s is empty \n',nf{k2});
        end
    end
    if misCrop > 0
        fprintf('Number of uncropped Stims:  %d \n',misCrop);
    else
        if nanStim > 0
            fprintf('Crop done completely with %d empty stimulation \n',nanStim);
        else
            disp('Crop done completely');
        end
    end
else
    disp('Crop done completely');
end

fprintf('\n\nNumber of Stims notes found: %d\n',realStims);
fprintf('Number of stimulations detected: %d\n', stimTxt);
fprintf('Number of stimulations cropped: %d\n\n',realCrps);
txtf = dir(fullfile(DirOut,'*.txt')); % delete all .txt file after cropping
txtf = {txtf.name};
for k3 = 1:numel(txtf)
    if isempty(strfind(txtf{k3},'CropsNoRecordings_log'))
        delete(fullfile(DirOut, txtf{k3}))
    end
end
if matDelCount > 0
    logTxtFile = fullfile(DirOut, [matFile,'CropsNoRecordings_log.txt']);
    fprintf('\n\n stim channels not recorded, Crop is deleted: \nn');
    try
        fid = fopen(logTxtFile,'a');
        for i = 1:matDelCount
            fprintf('%s \n', matDeleted{i});
            fprintf(fid,'%s\n',matDeleted{i});
            if exist(fullfile(DirOut, strcat(matDeleted{i},'.mat')), 'file') == 2
                delete(fullfile(DirOut, strcat(matDeleted{i},'.mat')))
                delete(fullfile(DirOut, strcat(matDeleted{i},'.dat')))
            end
        end
        fclose(fid);
    catch
        disp('Log not created.')
    end
end
set_final_status('OK');
disp('Done');

