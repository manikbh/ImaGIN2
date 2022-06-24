function Dnew = ImaGIN_BipolarMontage(S)
% Perform a longitudinal bipolar montage

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
% Authors: Olivier David

try
    Filename = S.Fname;
catch
    Filename = spm_select(Inf, '\.mat$', 'Select data file');
end

try
    FileOut=S.FileOut;
catch
    S.FileOut = Filename;
end

try
    Bipole=S.Bipole;
catch
    Bipole = [];
end

try
    SaveChannels = S.Save;
catch
    SaveChannels = 0;
end



for i0=1:size(Filename,1)
    try
        bipole = load(Bipole);
    catch
        bipole = Bipole;
    end
    T=deblank(Filename(i0,:));
    D=spm_eeg_load(T);
    
    Sensors=sensors(D,'EEG');
    Name=Sensors.label;
    Position=Sensors.elecpos;
    BadChannelsMono=badchannels(D);
    idxNaN = [] % find(isnan(Position(:,1))); We now assume channels without coordinates are not bad (CRFs older than 2009 do not have T1pre coordinates)
    if ~isempty(idxNaN)
        BadChannelsMono = [BadChannelsMono(:);idxNaN(:)];        
    end
    BadChannelsMono = sort(unique(BadChannelsMono));
    if isempty(bipole)
        Cpos2full=[];
        Cnamesfull={};
        bipole=[];       
        for i1=1:numel(Name)
            matched_1 = regexp(Name{i1},'^(?<label>.*?)(?<index>[0-9]+)$','once','names');
            if isempty(matched_1)
                continue
            end            
            label_name_1 = matched_1.label;
            label_number_1 = matched_1.index;
            matched_channels = 0;
            for i2=1:numel(Name)        
                matched_2 = regexp(Name{i2},'^(?<label>.*?)(?<index>[0-9]+)$','once','names');
                if isempty(matched_2)
                    continue
                end                 
                label_name_2 = matched_2.label;
                label_number_2 = matched_2.index;       
                if strcmp(label_name_1,label_name_2) & str2num(label_number_1)+1 == str2num(label_number_2) 
                    matched_channels = matched_channels+1;
                    Cpos2full(:,end+1)= mean(Position([i1 i2],:))';
                    Cnamesfull{1,end+1}=[Name{i1} '-' Name{i2}];
                    bipole =[bipole [i2;i1]];
                end        
            end
            if matched_channels > 1
                error('Multiple monopolar channels with the same label found.')
            end
        end
        if strcmp(Name{end},'ecg')
            Cpos2full(:,end+1) = mean(Position(end,:),1)';
            Cnamesfull{1,end+1} = Name{end};
            bipole =[bipole [length(Name);length(Name)]];
        end       
    else
        bipole=bipole';
        Cpos2full=zeros(3,size(bipole,2));
        Cnamesfull={};
        for i1=1:size(bipole,2)
            Cpos2full(:,i1)=mean(Position(bipole(:,i1),:),1)';
            Cnamesfull{1,end+1}=[Name{bipole(1,i1)} '-' Name{bipole(2,i1)}];
        end
    end

    
    Sensors.chanpos=Cpos2full';
    Sensors.elecpos=Cpos2full';
    Sensors.label=Cnamesfull;
    Sensors.chanunit=Sensors.chanunit(bipole(1,:));
    Sensors.chantype=Sensors.chantype(bipole(1,:));
    
    %New bad channels
    BadChannelsBip=[];
    for i1=1:length(BadChannelsMono)
        [a1,a2]=find(bipole==BadChannelsMono(i1));
        BadChannelsBip=[BadChannelsBip a2'];
    end
    BadChannelsBip=unique(BadChannelsBip);    
    
    
    
    Data=D(bipole(1,:),:,:)-D(bipole(2,:),:,:);
    
    %in case monopolar is kept on some electrodes
    Data(find(bipole(1,:)==bipole(2,:)),:,:)=D(bipole(1,find(bipole(1,:)==bipole(2,:))),:,:);
    
    if strcmp(Name{end},'ecg')
        Data(end,:)=D(bipole(1,end),:,:);
    end
    
    %Save as a newfile
    try
        newname=FileOut;
        Dnew=clone(D,newname, [size(Data,1) D.nsamples D.ntrials]);
    catch
        newname=['b' D.fname];
        Dnew=clone(D,newname, [size(Data,1) D.nsamples D.ntrials]);
    end
    Dnew(:,:,:)=Data;
    Dnew=sensors(Dnew,'EEG',Sensors);
    Dnew=chanlabels(Dnew,1:size(Data,1),Sensors.label);
    Dnew=chantype(Dnew,1:size(Data,1),chantype(D,bipole(1,:)));
    if ~isempty(BadChannelsBip)
        Dnew = badchannels(Dnew,BadChannelsBip,1); % add badchannel index in meeg object
    end
    
    save(Dnew);
    
    if SaveChannels
        [dir,name,~]=fileparts(newname);
        fid = fopen(fullfile(dir, ['recordings_' name '.txt']),'w');
        recording_output = (Dnew.sensors('eeg').label);
        recording_outputGood = setdiff(recording_output, recording_output(BadChannelsBip));
        fprintf(fid,'%s\n',recording_outputGood{:});
        fclose(fid);
    end
    
    % Save the list of channels in the output file
    ImaGIN_save_log(fullfile(Dnew), 'Bipolar montage:', chanlabels(Dnew));
end

