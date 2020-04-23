function prepare_ImaGIN_SpikesDetection(path_in,stims_nb,varargin)
% Dummy version of prepare_ImaGIN_SpikesDetection to test variable ammount
% of inputs passed from the pipeline
    fprintf('path_in is: %s \n', path_in)
    fprintf('stims_nb is: %s \n', stims_nb)
    for ii=1:numel(varargin)
        fprintf('varargin %d is: %s \n',ii,varargin{ii})
    end
    fprintf('Done\n')
end