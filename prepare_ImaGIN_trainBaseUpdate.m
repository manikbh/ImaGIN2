function prepare_ImaGIN_trainBaseUpdate(inDir)
% FileIn: path linking to the Cropped MEEG object
clear S;
S = fullfile(inDir);
ImaGIN_trainBaseUpdate(S);  
close all;
end

