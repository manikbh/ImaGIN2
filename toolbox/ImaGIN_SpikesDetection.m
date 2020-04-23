function ImaGIN_SpikesDetection(S)
% This function extracts and concatenates the baselines from cropped, qualitatively controlled, stimulation files.
% It runs Delphos' spike detection on extracted baselines after bipolar montage computation.
% It then computes corrected spike rates for each bipolar recording channel.

    % Vars parsing
    stims_files = S.stims_files;
    path_out = S.path_out;

    % Extract baselines: get the pre-stimulation period for each stimulation and concatenate them for each channel.
    concat_baselines = []; % Store timeseries
    concat_zeroed = []; % Store booleans later used as logic gates on timeseries
    for ff=1:numel(stims_files)
        D = spm_eeg_load(stims_files{ff});  
        time_axis = time(D);
        onset = timeonset(D);    
        negative_time_idxes = find(time_axis > onset+1 & time_axis < -1); % Arbitrary time margins of 1 s around prestimulus period    
        % Extract baseline of the stimulation
        file_baselines = D(:,negative_time_idxes,1); 
        concat_baselines = [concat_baselines file_baselines]; 
        % Generate the corresponding boolean matrix to identify good (1) versus bad (0) samples
        file_zeroed = ones(size(file_baselines));
        all_channels_idxes = indchantype(D,'ALL');
        good_eeg_channels_idxes = indchantype(D,'EEG','GOOD');          
        channels_to_zeroed = setdiff(all_channels_idxes,good_eeg_channels_idxes); 
        file_zeroed(channels_to_zeroed,:,1) = zeros(numel(channels_to_zeroed),numel(negative_time_idxes));
        concat_zeroed = [concat_zeroed file_zeroed];
    end

    % Generate the baseline object.
    baseline_file = fullfile(path_out,'baseline');
    baseline_obj = clone(D,baseline_file,[size(concat_baselines,1) size(concat_baselines,2) 1],3);
    baseline_obj(:) = concat_baselines(:);
    baseline_obj = events(baseline_obj,[],{''});
    baseline_obj = timeonset(baseline_obj,0);
    save(baseline_obj);
    
    % Compute bipolar montage on baseline object to get the bipolar baseline object
    S.Fname = baseline_file;
    bp_baseline_obj = ImaGIN_BipolarMontage(S);  
    save(bp_baseline_obj);
    
    % Estimate bipolar logic gates from monopolar boolean matrix   
    bp_concat_zeroed = zeros(size(bp_baseline_obj));
    for cc=1:nchannels(bp_baseline_obj)
        bp_chan = chanlabels(bp_baseline_obj,cc);
        mono_chans = regexp(bp_chan,'-','split'); % split(bp_chan,'-')
        chan1_idx = indchannel(D,mono_chans{1}{1}); % mono_chans{1}
        chan2_idx = indchannel(D,mono_chans{1}{2}); % mono_chans{2}
        both_good = (concat_zeroed(chan1_idx,:,1) == 1 & concat_zeroed(chan2_idx,:,1) == 1); 
        bp_concat_zeroed(cc,find(both_good)) = ones(1,numel(find(both_good))); 
    end

    % Multiply the bipolar baselines with the bipolar logic gates and create the gated bipolar baseline object
    gated_bp_timeseries = [];
    for cc=1:nchannels(bp_baseline_obj) % This loop exists because of memory issues.
        gated_bp_timeseries(cc,:) = bp_baseline_obj(cc,:).*bp_concat_zeroed(cc,:);        
    end   
    gated_bp_baseline_obj = clone(bp_baseline_obj,fullfile(path_out,'gated_bp_baseline')); 
    gated_bp_baseline_obj(:,:) = gated_bp_timeseries;
    save(gated_bp_baseline_obj);

    % Run Delphos' spike detection on gated baselines
    fs = gated_bp_baseline_obj.fsample; 
    results = Delphos_detector(gated_bp_baseline_obj,chanlabels(gated_bp_baseline_obj),'SEEG',fs,{'Spk'},[],[],40,[]);

    % Calculate and store corrected spike rates in a table
    total_samples = nsamples(gated_bp_baseline_obj);
    total_times = [];    
    for cc=1:nchannels(gated_bp_baseline_obj)
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
    [~, bname,~] = fileparts(gated_bp_baseline_obj.fnamedat);
    spkFileName = ['SPK_' bname '.mat'];
    spkFile.spikeRate = rates_table;
    spkFile.markers   = spks_table;
    spkFile.baseline  = ['Total duration of ' num2str((nsamples(gated_bp_baseline_obj)/fsample(gated_bp_baseline_obj))/60)  ' minutes'];
    spkFile.comment   = 'Spike rate is number of spikes per minute whithin a bipolar channel';
    save(fullfile(path_out,spkFileName), 'spkFile');
    
    % Store rate table in an individual files for easy extraction
    rate_table_file = ['spike_rate_' bname '.mat'];
    spike_rates = [channels,num2cell(rates)];
    save(fullfile(path_out,rate_table_file), 'spike_rates');
    
    % Ending step
    set_final_status('OK');
    disp('Done');
end
