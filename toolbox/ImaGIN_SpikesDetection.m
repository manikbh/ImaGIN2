function ImaGIN_SpikesDetection(S)
% This function extracts and concatenates the baselines from cropped, qualitatively controlled, stimulation files.
% It runs Delphos' spike detection on extracted baselines after bipolar montage computation.
% It then computes corrected spike rates for each bipolar recording channel.
% /!\ Cropped files CAN have a variable number of channels and thus the monopolar <-> bipolar mapping has to be established for each cropped file
% and taken into account during baseline concatenation for consistency. /!\

    %% Vars parsing
    stims_files = S.stims_files;
    path_out = S.path_out; 
    channels = struct;
    
    %% Extract data from monopolar stimulation files  
    for ff=1:numel(stims_files)
        % Establish the mappping for the stimulation
        stim_mapping = monopolar_to_bipolar_mapping(stims_files{ff},path_out);
        % Load the stimulation for data extraction
        D = spm_eeg_load(stims_files{ff});
        % Get the samping frequency
        sampling(ff) = fsample(D);
        % Get the coordinates
        sensors_info = sensors(D,'EEG');
        % Get the list of bad monopolar channels
        all_channels_idxes = indchantype(D,'ALL');
        good_eeg_channels_idxes = indchantype(D,'EEG','GOOD');          
        mono_bad_channel_idxes = setdiff(all_channels_idxes,good_eeg_channels_idxes);
        % Use the mapping to identify monopolar channels later coupled with bad channels
        bp_bad_channels_idxes = [];
        for cc=1:numel(mono_bad_channel_idxes)           
            bad_channel_1_map = find(mono_bad_channel_idxes(cc) == [stim_mapping.mono_idx1]);
            bad_channel_2_map = find(mono_bad_channel_idxes(cc) == [stim_mapping.mono_idx2]);            
            circular_bp_bad_channels_idxes = unique([stim_mapping(bad_channel_1_map).mono_idx1,...
            stim_mapping(bad_channel_1_map).mono_idx2,...
            stim_mapping(bad_channel_2_map).mono_idx1,...
            stim_mapping(bad_channel_2_map).mono_idx2]);       
            bp_bad_channels_idxes = unique([bp_bad_channels_idxes circular_bp_bad_channels_idxes]);           
        end
        % Concatenate all identified bad channels
        all_bad_channels_idxes = unique([mono_bad_channel_idxes bp_bad_channels_idxes]);
        % Extract baseline of the stimulation and zero all bad channels
        time_axis = time(D);
        onset = timeonset(D);        
        negative_time_idxes = find(time_axis > onset+1 & time_axis < -1); % Arbitrary time margins of 1 s around prestimulus period              
        file_baselines = D(:,negative_time_idxes,1);        
        file_baselines(all_bad_channels_idxes,:,1) = zeros(numel(all_bad_channels_idxes),numel(negative_time_idxes));  
        % Generate the corresponding boolean matrix to identify good (1) versus bad (0) samples
        file_zeroed = ones(size(file_baselines));
        file_zeroed(all_bad_channels_idxes,:,1) = zeros(numel(all_bad_channels_idxes),numel(negative_time_idxes));
        % Update the channels struct 
        eeg_chans_idxes = indchantype(D,'EEG');
        for cc=1:numel(eeg_chans_idxes)
            channel_name = chanlabels(D,eeg_chans_idxes(cc));
            channels.(channel_name{:}).timeserie{ff,:} = file_baselines(eeg_chans_idxes(cc),:,:);  
            channels.(channel_name{:}).zeroed{ff,:} = file_zeroed(eeg_chans_idxes(cc),:,:);            
            channels.(channel_name{:}).coordinates{ff,:} = sensors_info.elecpos(find(strcmp(channel_name,sensors_info.label)),:);            
        end
        % Tidy the channels struct with empty inputs for consistency
        temp_montage = fieldnames(channels);
        for nn=1:numel(temp_montage)
            if sum(strcmp(temp_montage{nn},chanlabels(D))) == 0
                channels.(temp_montage{nn}).timeserie{ff} = [];
                channels.(temp_montage{nn}).zeroed{ff} = [];
                channels.(temp_montage{nn}).coordinates{ff} = [];
            end
        end        
    end
    
    %% Test for extracted baseline sizes and fill-in missing channels
    global_montage = fieldnames(channels);
    for ff=1:numel(stims_files)        
        for nn=1:numel(global_montage)
            series_length(nn) = numel(channels.(global_montage{nn}).timeserie{ff,:});            
        end
        empty_series = find(series_length == 0);        
        nonempty_series = find(series_length);        
        valid_length = series_length(nonempty_series);        
        if all(valid_length == valid_length(1))            
            for ee=1:numel(empty_series)
                channels.(global_montage{empty_series(ee)}).timeserie{ff,:} = zeros(1,valid_length(1));
                channels.(global_montage{empty_series(ee)}).zeroed{ff,:} = zeros(1,valid_length(1));
                channels.(global_montage{empty_series(ee)}).coordinates{ff,:} = NaN(1,3);
            end  
        else
            error('Inconsistent length of channels.');
        end    
    end   
    
    %% Create the new monopolar baseline object
    for nn=1:numel(global_montage)
        channel{nn} = global_montage{nn};
        timeseries(nn,:) =  [channels.(global_montage{nn}).timeserie{:}];
        zeroed(nn,:) =  [channels.(global_montage{nn}).zeroed{:}];
        coordinates(nn,:) = channels.(global_montage{nn}).coordinates{1,:};
    end   
    % Create the new sensor
    new_sensors_info.balance.current = 'none';
    new_sensors_info.chanpos = coordinates;
    new_sensors_info.chantype = cell(1,numel(global_montage)); new_sensors_info.chantype(:) = {'eeg'};
    new_sensors_info.chanunit = cell(1,numel(global_montage)); new_sensors_info.chanunit(:) = {'V'};
    new_sensors_info.elecpos = coordinates;
    new_sensors_info.label = channel;
    new_sensors_info.type = 'ctf';
    new_sensors_info.unit = 'mm';
    % Create the new files 
    file_header = fullfile(path_out,'monopolar_baseline.mat');
    file_data = fullfile(path_out,'monopolar_baseline.dat');    
    if exist(file_header, 'file') == 2
     delete(file_header);      
    end
    if exist(file_data, 'file') == 2
     delete(file_data);          
    end   
    % Create the new object    
    D = meeg(size(timeseries,1),size(timeseries,2),1);
    D = fname(D,'monopolar_baseline');
    D = path(D,path_out);
    if all(sampling == sampling(1))
        D = fsample(D,sampling(1));
    else
        error('Inconsistent sampling frequency across cropped files.');
    end
    D = sensors(D,'EEG',new_sensors_info)    
    D = blank(D) ;
    D = link(D,file_data);   
    D(:,:,1) = timeseries;   
    save(D);    
    
    %% Compute bipolar montage on baseline object to get the bipolar baseline object
    Sbp.Fname = file_header;
    Sbp.FileOut = fullfile(path_out,'bipolar_baseline');
    bp_baseline_obj = ImaGIN_BipolarMontage(Sbp);     
    global_mapping = monopolar_to_bipolar_mapping(Sbp.Fname,path_out);

    % Estimate bipolar boolean matrix from monopolar boolean matrix   
    bp_concat_zeroed = zeros(size(bp_baseline_obj));
    for cc=1:numel(global_mapping)  
        both_good = (zeroed(global_mapping(cc).mono_idx1,:,1) & zeroed(global_mapping(cc).mono_idx2,:,1)); 
        bp_concat_zeroed(cc,find(both_good)) = ones(1,numel(find(both_good))); 
    end
    
    % Run Delphos' spike detection on baselines
    fs = bp_baseline_obj.fsample; 
    results = Delphos_detector(bp_baseline_obj,chanlabels(bp_baseline_obj),'SEEG',fs,{'Spk'},[],[],40,[]);

    % Calculate and store corrected spike rates in a table
    total_samples = nsamples(bp_baseline_obj);
    total_times = [];    
    for cc=1:nchannels(bp_baseline_obj)
        zeroed_samples = numel(find(~bp_concat_zeroed(cc,:,1)));      
        total_times(cc) = (total_samples - zeroed_samples) /fs/60; % Total recording time in minutes
    end
    chans = [results.labels]';
    rates = round(results.n_Spk./total_times',1); % Round to nearest decimal
    rates_table  = table(chans,rates);

    % Store spike features in a table
    spks_positions = [results.markers.position]';
    spks_type = {results.markers.label}';
    spks_channels = {results.markers.channels}';
    spks_values = {results.markers.value}';
    spks_table = table(spks_channels,spks_type,spks_positions,spks_values);

    % Store all tables in a .mat
    [~,bname,~] = fileparts(bp_baseline_obj.fname);
    spkFileName = ['Results_' bname '.mat'];
    spkFile.spikeRate = rates_table;
    spkFile.markers   = spks_table;
    spkFile.baseline  = ['Total duration of ' num2str((nsamples(bp_baseline_obj)/fsample(bp_baseline_obj))/60)  ' minutes'];
    spkFile.comment   = 'Spike rate is number of spikes per minute whithin a bipolar channel';
    save(fullfile(path_out,spkFileName), 'spkFile');
    
    % Store rate table in an individual files for easy extraction
    rate_table_file = 'SPK_bBaseline_rates.mat';
    spike_rates = [chans,num2cell(rates)];
    save(fullfile(path_out,rate_table_file), 'spike_rates');
    
    % Ending step
    set_final_status('OK');
    disp('Done');
end

function channels_map = monopolar_to_bipolar_mapping(monopolar_file,path_out)
    mono_stimulation_obj = spm_eeg_load(monopolar_file);
    monopolar_sensors = sensors(mono_stimulation_obj,'EEG');
    [mono_path,mono_name,mono_ext] = fileparts(monopolar_file);
    Sbp.Fname = monopolar_file;
    Sbp.FileOut = fullfile(path_out,['b' mono_name]);
    bp_stimulation_obj = ImaGIN_BipolarMontage(Sbp);     
    for cc=1:nchannels(bp_stimulation_obj)
        bipolar_chan = chanlabels(bp_stimulation_obj,cc);   
        mono_chans = regexp(bipolar_chan,'-','split');  
        channels_map(cc).bp_chan = bipolar_chan;
        channels_map(cc).bp_idx = cc;
        channels_map(cc).mono_chan1 = mono_chans{1}{1};
        channels_map(cc).mono_chan2 = mono_chans{1}{2};        
        % Get the monopolar channels idxes.
        chan1_idx = find(strcmp(monopolar_sensors.label,mono_chans{1}{1})); 
        chan2_idx = find(strcmp(monopolar_sensors.label,mono_chans{1}{2}));        
        % Test and correct for duplicated labels in the monopolar object
        if numel(chan1_idx) == 0
            error('Could not find the monopolar channel from the bipolar montage.');
        elseif numel(chan1_idx) == 1 % If 1 match, as expected, do nothing
            %pass        
        else
            error('2 or more channels with the same label')        
        end
        % Same logic as for chan1_idx
        if numel(chan2_idx) == 0
            error('Could not find the monopolar channel from the bipolar montage.');
        elseif numel(chan2_idx) == 1
            %pass       
        else
            error('2 or more channels with the same label')        
        end
        channels_map(cc).mono_idx1 = chan1_idx;
        channels_map(cc).mono_idx2 = chan2_idx;        
    end    
    delete([Sbp.FileOut '.mat']);delete([Sbp.FileOut '.dat']);delete([Sbp.FileOut '_log.txt']); 
end
