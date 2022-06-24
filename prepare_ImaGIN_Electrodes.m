function prepare_ImaGIN_Electrodes(source, polarity, FileIn, implantationFile, FileOut, FileTxtOut)
 
%Add electrode position
%source : 'classic' (a file with the electrodes names and another one with
%the MNI positions), implantationFile is the source directory of the 2
%files
%      or 'intranat' (a file -implantationFile- with electrodes names, position, atlas,...)

[Root,file,~] = fileparts(FileIn);

clear S
S.Fname = fullfile(Root,[spm_str_manip(file,'r') '.mat']);
if strcmp(polarity,'monopolar')
    if strcmp(source,'classic')
        S.filenamePos = fullfile(implantationFile,'Electrodes_Pos_MNI.txt');
        S.filenameName = fullfile(implantationFile,'Electrodes_Name.txt');
    elseif strcmp(source, 'intranat')
        S.filenameName = implantationFile;
    end
else
    error('No other polarity option is implemented here!')
end

% %for files generated with Intranat
% if strcmp(source, 'intranat')
%     data=extractIntranatFile(implantationFile);
%     for i1=1:size(data,2)
%         for i2=1:size(data{i1},1)
%             if strcmp(data{i1}{i2},'MNI')
%                 indi=i1;
%                 indj=i2;
%             end
%         end
%     end
%     positiontmp = data{indi};
%     nametmp = data{1};
%     Position = [];
%     Name = cell(1,1);
%     if strcmp(polarity,'monopolar')
%         for i=(indj+1):length(nametmp)
%             if isempty(strfind(nametmp{i},'-'))
%                 
%                 iLastLetter = find(~ismember(nametmp{i}, '0123456789'), 1, 'last');
%                 if isempty(iLastLetter)
%                     iLastLetter = 0;
%                 end
%                 
%                 if iLastLetter > 0
%                     name     = [lower(nametmp{i}(1:iLastLetter)) nametmp{i}(iLastLetter+1:end)];
%                     Name{1}  = cat(1,Name{1},{name});
%                     Position = cat(1,Position,str2num(positiontmp{i}));
%                 end
%             end
%         end
%     else
%         if strcmp(polarity, 'bipolar')
%             % TODO: Throw an exception here, don't execute anything!!!
%             error('Bipolar option not implemented')
%         end
%     end
% S.Position = Position;
% S.Name = Name;
% end
% 

S.FileTxtOut = FileTxtOut;
S.FileOut = FileOut;
fprintf('prepare_ImaGIN_Electrodes: S.FileTxtOut: %s \n', S.FileTxtOut)
fprintf('prepare_ImaGIN_Electrodes: S.FileOut: %s \n', S.FileOut)

ImaGIN_Electrode(S);

set_final_status('OK')

close all

end

% % It gets the CSV file, and reads its content to Matlab structure 'data'
% function data = extractIntranatFile(filename)
% delimiter = '\t';
% formatSpec = '%s%s%s%s%s%s%s%s%s%s%s%s%s%s%[^\n\r]';
% 
% fileID = fopen(filename,'r');
% data = textscan(fileID, formatSpec, 'Delimiter', delimiter,  'ReturnOnError', false);
% fclose(fileID);
% end
