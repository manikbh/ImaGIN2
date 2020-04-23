function prepare_ImaGIN_SpikesDetection(path_in,stims_nb,varargin)
% Concatenates baselines from individual stimulation files and runs
% delphos' spike detection to compute spike rates for each bipolar channel.

    % Input parameter parsing   
    p = inputParser; 
    addRequired(p,'path_in',@ischar);
    addRequired(p,'stims_nb',@ischar); 
    parse(p,path_in,stims_nb); 
    
    % Prepare vars for ImaGIN_SpikesDetection 
    path_in = p.Results.path_in;
    stims_nb = str2num(p.Results.stims_nb);
    S.path_out = varargin{end};
    S.stims_files = {};
    for ii=1:stims_nb
        S.stims_files{ii} = [path_in filesep varargin{ii}];
    end
    
    % Call to ImaGIN_SpikesDetection 
    fprintf('prepare_ImaGIN_SpikesDetection: Directory: %s \n', path_in)
    ImaGIN_SpikesDetection(S);
end