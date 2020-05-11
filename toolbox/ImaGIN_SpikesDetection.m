function ImaGIN_SpikesDetection(S)
% This function extracts and concatenates the baselines from cropped, qualitatively controlled, stimulation files.
% It runs Delphos' spike detection on extracted baselines after bipolar montage computation.
% It then computes corrected spike rates for each bipolar recording channel.
% Manipulating the data of hundreds of stimulation files brings performance (RAM) issues...
% ...which are addressed by zeroing the baselines during the extraction loop, at the expense of readability.

    % Vars parsing
    stims_files = S.stims_files;
    path_out = S.path_out;
    
    % Establish the channel mapping from monopolar to bipolar using a single stimulation file    
    mono_stimulation_file = fullfile(path_out,'single_stimulation');
    mono_stimulation_obj = clone(spm_eeg_load(stims_files{1}),mono_stimulation_file); 
    Sbp.Fname = mono_stimulation_file;
    bp_stimulation_obj = ImaGIN_BipolarMontage(Sbp);   
    clear channels_map;
    for cc=1:nchannels(bp_stimulation_obj)
        bipolar_chan = chanlabels(bp_stimulation_obj,cc);
        mono_chans = regexp(bipolar_chan,'-','split');  
        channels_map(cc).bp_chan = bipolar_chan;
        channels_map(cc).bp_idx = cc;
        channels_map(cc).mono_chan1 = mono_chans{1}{1};
        channels_map(cc).mono_chan2 = mono_chans{1}{2};        
        
        % Test if the channel label exists in the monopolar object, if not, add a '0' and try again.
        chan1_idx = indchannel(mono_stimulation_obj,mono_chans{1}{1});
        if numel(chan1_idx) == 0
            old_name = mono_chans{1}{1};
            match = regexp(old_name,'^(?<letters>[a-zA-Z]+)(?<digits>[0-9]+)$','once','names');
            new_name = [match.letters '0' match.digits];
            channels_map(cc).mono_chan1 = new_name;
            chan1_idx = indchannel(mono_stimulation_obj,new_name);
        end       
        chan2_idx = indchannel(mono_stimulation_obj,mono_chans{1}{2});
        if numel(chan2_idx) == 0
            old_name = mono_chans{1}{2};
            match = regexp(old_name,'^(?<letters>[a-zA-Z]+)(?<digits>[0-9]+)$','once','names');
            new_name = [match.letters '0' match.digits];
            channels_map(cc).mono_chan2 = new_name;
            chan2_idx = indchannel(mono_stimulation_obj,new_name);
        end  
        
        % Test and correct for duplicated labels in the monopolar object
        if numel(chan1_idx) == 0 % Could not find the monopolar label
            error('Could not find the monopolar channel from the bipolar montage.');
        elseif numel(chan1_idx) == 1 % If 1 match, as expected, do nothing
            %pass
        elseif numel(chan1_idx) == 2 % If 2 matches, find and eliminate the channel without coordinates AND flagged as badchannel (scalp measurement)
            monopolar_sensors = sensors(mono_stimulation_obj,'EEG');            
            chan1_matching = [sum(badchannels(mono_stimulation_obj) == chan1_idx(1)) == 1 & sum(isnan(monopolar_sensors.elecpos(chan1_idx(1),:))) == 3,...
                sum(badchannels(mono_stimulation_obj) == chan1_idx(2)) == 1 & sum(isnan(monopolar_sensors.elecpos(chan1_idx(2),:))) == 3]; % [0,1] or [1,0]
            chan1_idx = chan1_idx(~chan1_matching(:));
            clear monopolar_sensors;
        else
            error('3 or more channels with the same label')        
        end
        % Same logic as for chan1_idx
        if numel(chan2_idx) == 0
            error('Could not find the monopolar channel from the bipolar montage.');
        elseif numel(chan2_idx) == 1
            %pass
        elseif numel(chan2_idx) == 2
            monopolar_sensors = sensors(mono_stimulation_obj,'EEG');            
            chan2_matching = [sum(badchannels(mono_stimulation_obj) == chan2_idx(1)) == 1 & sum(isnan(monopolar_sensors.elecpos(chan2_idx(1),:))) == 3,...
                sum(badchannels(mono_stimulation_obj) == chan2_idx(2)) == 1 & sum(isnan(monopolar_sensors.elecpos(chan2_idx(2),:))) == 3];            
            chan2_idx = chan2_idx(~chan2_matching(:));
            clear monopolar_sensors;
        else
            error('3 or more channels with the same label')        
        end     
        
        channels_map(cc).mono_idx1 = chan1_idx;
        channels_map(cc).mono_idx2 = chan2_idx;
    end
    clear mono_stimulation_obj bp_stimulation_obj; 
    
    % Extract baselines: get the pre-stimulation period for each stimulation and concatenate them for each channel.
    concat_baselines = []; % Store timeseries
    concat_zeroed = []; % Store booleans later used for the correction of spike rates
    for ff=1:numel(stims_files)
        clear D;
        D = spm_eeg_load(stims_files{ff});                
        % Get the list of monopolar bad channels
        all_channels_idxes = indchantype(D,'ALL');
        good_eeg_channels_idxes = indchantype(D,'EEG','GOOD');          
        mono_bad_channel_idxes = setdiff(all_channels_idxes,good_eeg_channels_idxes);
        % Use the monopolar-bipolar channel mapping to identify monopolar channels later coupled with bad channels
        bp_bad_channels_idxes = [];
        for cc=1:numel(mono_bad_channel_idxes)           
            bad_channel_1_map = find(mono_bad_channel_idxes(cc) == [channels_map.mono_idx1]);
            bad_channel_2_map = find(mono_bad_channel_idxes(cc) == [channels_map.mono_idx2]);            
            circular_bp_bad_channels_idxes = unique([channels_map(bad_channel_1_map).mono_idx1,...
            channels_map(bad_channel_1_map).mono_idx2,...
            channels_map(bad_channel_2_map).mono_idx1,...
            channels_map(bad_channel_2_map).mono_idx2]);       
            bp_bad_channels_idxes = unique([bp_bad_channels_idxes circular_bp_bad_channels_idxes]);           
        end
        % Concatenate all identified bad channels
        all_bad_channels_idxes = unique([mono_bad_channel_idxes bp_bad_channels_idxes]);
        % Extract baseline of the stimulation and zeroed all bad channels
        time_axis = time(D);
        onset = timeonset(D);        
        negative_time_idxes = find(time_axis > onset+1 & time_axis < -1); % Arbitrary time margins of 1 s around prestimulus period    
        file_baselines = D(:,negative_time_idxes,1);        
        file_baselines(all_bad_channels_idxes,:,1) = zeros(numel(all_bad_channels_idxes),numel(negative_time_idxes));       
        concat_baselines = [concat_baselines file_baselines];        
        % Generate the corresponding boolean matrix to identify good (1) versus bad (0) samples
        file_zeroed = ones(size(file_baselines));
        file_zeroed(all_bad_channels_idxes,:,1) = zeros(numel(all_bad_channels_idxes),numel(negative_time_idxes));
        concat_zeroed = [concat_zeroed file_zeroed];
    end

    % Generate the baseline object.
    baseline_file = fullfile(path_out,'Baseline');
    baseline_obj = clone(D,baseline_file,[size(concat_baselines,1) size(concat_baselines,2) 1],3);
    baseline_obj(:) = concat_baselines(:);
    baseline_obj = events(baseline_obj,[],{''});
    baseline_obj = timeonset(baseline_obj,0);
    clear baseline_obj;
    
    % Compute bipolar montage on baseline object to get the bipolar baseline object
    S.Fname = baseline_file;
    bp_baseline_obj = ImaGIN_BipolarMontage(S);  
    
    % Estimate bipolar boolean matrix from monopolar boolean matrix   
    bp_concat_zeroed = zeros(size(bp_baseline_obj));
    for cc=1:numel(channels_map)  
        both_good = (concat_zeroed(channels_map(cc).mono_idx1,:,1) == 1 & concat_zeroed(channels_map(cc).mono_idx2,:,1) == 1); 
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
    channels = [results.labels]';
    rates = round(results.n_Spk./total_times',1); % Round to nearest decimal
    rates_table  = table(channels,rates);

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
    spike_rates = [channels,num2cell(rates)];
    save(fullfile(path_out,rate_table_file), 'spike_rates');
    
    % Ending step
    set_final_status('OK');
    disp('Done');
end
