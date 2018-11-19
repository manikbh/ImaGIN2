function prepare_ImaGIN_Electrodes(source, polarity, FileIn, implantationFile, FileOut, FileTxtOut)
 
%Add electrode position
%source : 'classic' (a file with the electrodes names and another one with
%the MNI positions), implantationFile is the source directory of the 2
%files
%      or 'intranat' (a file -implantationFile- with electrodes names, position, atlas,...)

[Root,file,~] = fileparts(FileIn);

clear S
S.Fname = fullfile(Root,[spm_str_manip(file,'r') '.mat']);

if strcmp(source,'classic')
    S.filenamePos = fullfile(implantationFile,'Electrodes_Pos_MNI.txt');
    S.filenameName = fullfile(implantationFile,'Electrodes_Name.txt');
end

%for files generated with Intranat
if strcmp(source, 'intranat')
    % data will contain the content of the implantationFile (presumably .CSV)
    data=extractIntranatFile(implantationFile);
    
    %identification of the line with the MNI positions
    for i1=1:size(data,2)
        for i2=1:size(data{i1},1)
            if strcmp(data{i1}{i2},'MNI')
                % getting indeces of the 'MNI' column
                indi=i1;
                indj=i2;
            end
        end
    end
    positiontmp = data{indi};
    nametmp = data{1};
    Position = [];
    Name = cell(2,1);
    if strcmp(polarity,'monopolar')
        monop=[];
        for i=(indj+1):length(nametmp)
            if isempty(strfind(nametmp{i},'-'))
                ind1 = regexp(lower(nametmp{i}),'[a-z\'']');
                ind2 = regexp(lower(nametmp{i}),'[0-9]');
                select_ind=[];
                for i0=1:length(ind2)
                    if ind2(i0)>ind1(length(ind1))
                        select_ind=cat(1,select_ind,ind2(i0));
                    end
                end
                % separate contact label and contact index
                num = num2str(str2num(nametmp{i}(select_ind)));
                num2 = nametmp{i}(select_ind);
                name = [lower(nametmp{i}(ind1(1):ind1(length(ind1)))) num];
                name2 = [lower(nametmp{i}(ind1(1):ind1(length(ind1)))) num2];
                Name{1} = cat(1,Name{1},{name});
                Name{2} = cat(1,Name{2},{name2});
                Position = cat(1,Position,str2num(positiontmp{i}));
            end
        end
        
    else
        if strcmp(polarity, 'bipolar')  % fixed even for the case where the electrode name includes numbers
        % TODO: Throw an exception here, don't execute anything!!!
            error('Bipolar option not implemented')
          %       
          %       iLastLetter1 = find(~ismember(nametmp{i}(1:sep-1),   '0123456789'), 1, 'last'); %VT
          %       iLastLetter2 = find(~ismember(nametmp{i}(sep+1:end), '0123456789'), 1, 'last');
          %       if isempty(iLastLetter1)
          %           iLastLetter1 = 0;
          %       end
          %       if isempty(iLastLetter2)
          %           iLastLetter2 = 0;
          %       end
          %       if iLastLetter1 > 0 && iLastLetter2 > 0
          %           name     = [lower(nametmp{i}(1:iLastLetter1)) nametmp{i}(iLastLetter1+1:sep-1) lower(nametmp{i}(1:iLastLetter1)) nametmp{i}(sep+iLastLetter2+1:end)];
          %           name2    = [lower(nametmp{i}(1:iLastLetter1)) nametmp{i}(iLastLetter1+1:sep-1) lower(nametmp{i}(1:iLastLetter1)) nametmp{i}(sep+iLastLetter2+1:end)];
          %           Name{1}  = cat(1,Name{1},{name});
          %           Name{2}  = cat(1,Name{2},{name2});
          %           Position = cat(1,Position,str2num(positiontmp{i}));
          %       end
          %       %{
          %       % VT
          %       %ind2 = regexp(nametmp{i},'\d');
          %       if ~isempty(ind1) && ~isempty(ind2)
          %           name = [lower(nametmp{i}(ind1 < sep)) num2str(str2num(nametmp{i}(ind2(ind2 < sep)))) lower(nametmp{i}(ind1 < sep)) num2str(str2num(nametmp{i}(ind2(ind2 > sep))))];
          %           name2 = [lower(nametmp{i}(ind1 < sep)) nametmp{i}(ind2(ind2 < sep)) lower(nametmp{i}(ind1 < sep)) nametmp{i}(ind2(ind2 > sep))];
          %           Name{1} = cat(1,Name{1},{name});
          %           Name{2} = cat(1,Name{2},{name2});
          %           Position = cat(1,Position,str2num(positiontmp{i}));
          %       end
          %       %}
          %   end
        end
        end
    end
S.Position = Position;
S.Name = Name;
end


S.FileTxtOut = FileTxtOut;
S.FileOut = FileOut;
fprintf('prepare_ImaGIN_Electrodes: S.FileTxtOut: %s \n', S.FileTxtOut)
fprintf('prepare_ImaGIN_Electrodes: S.FileOut: %s \n', S.FileOut)

D = ImaGIN_Electrode(S);

set_final_status('OK')

close all

end

% It gets the CSV file, and reads its content to Matlab structure 'data'
function data = extractIntranatFile(filename)
delimiter = '\t';
formatSpec = '%s%s%s%s%s%s%s%s%s%s%s%s%s%s%[^\n\r]';

fileID = fopen(filename,'r');
data = textscan(fileID, formatSpec, 'Delimiter', delimiter,  'ReturnOnError', false);
fclose(fileID);
end
