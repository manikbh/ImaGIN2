function prepare_ImaGIN_trainBaseUpdate(inDir)
% FileIn: path linking to the Cropped MEEG object
clear S;
S.DirIn = inDir;
ImaGIN_trainBaseUpdate(S);
   
f_id = fopen('/home/anthonyboyer/Documents/Debugs/20122019/ouput/log.txt','w');
fprintf(f_id,inDir);
fclose(f_id);

close all

end

