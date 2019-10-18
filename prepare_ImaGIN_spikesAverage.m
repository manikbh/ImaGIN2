function prepare_ImaGIN_spikesAverage(spkDirs)
clear S;
S.Dirs = fullfile(spkDirs); 
ImaGIN_spikesAverage(S);
set_final_status('OK')
close all
end
