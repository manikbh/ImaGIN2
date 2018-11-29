function prepare_ImaGIN_BadChannel(trainBase, toCheck, FileIn, FileOut)
   
% FileIn: path linking to the MEEG file to extract bad channel indices
clear S;
S.dataset = fullfile(FileIn); 
S.FileOut = FileOut;
S.trainBase = trainBase; % directory where a trained model resides
S.toCheck  = toCheck;
fprintf('prepare_ImaGIN_BadChannel: S.dataset: %s \n', S.dataset)
fprintf('prepare_ImaGIN_BadChannel: S.FileOut: %s \n', S.FileOut)
fprintf('prepare_ImaGIN_BadChannel: S.trainBase: %s \n', S.trainBase)
fprintf('prepare_ImaGIN_BadChannel: S.toCheck: %s \n', S.toCheck)

D = ImaGIN_BadChannel(S);

close all

end
