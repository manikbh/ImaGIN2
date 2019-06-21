
function newLabels = ImaGIN_channels_standard(allLabels, caseType, montageType)
    % Channel labels standardization: 'i''2i''1' becomes Ip02-Ip01
    if nargin < 3 || isempty(caseType), caseType = 'all_upper'; end    
    newLabels = {};        
    for i = 1: length(allLabels)
        Label = allLabels{i};
        % Remove spaces
        Label(Label == ' ') = [];
        % Switch case type
        switch (caseType)
            case 'no_change'
                % Nothing to change
            case 'all_upper'
                Label = upper(Label);
            case 'upper_except_p'
                    iP = find(Label == 'p');
                    Label = upper(Label);
                    if any(iP >= 2) && strcmpi(montageType,'monopolar')
                        Label(iP(end)) = 'p';
                    end
        end
        % Replacing ' with p
        Label = strrep(Label, '''', 'p');
        % Replacing , with p
        Label = strrep(Label, ',', 'p'); 
        
        [numb, idx] = regexp(Label,'\d*','Match');
        iLastLetter = find(~ismember(Label, '0123456789'), 1, 'last');
        chInd = Label(iLastLetter+1:end);
        if strcmpi(montageType,'bipolar')
            if numel(numb) == 2 && strcmp(chInd, numb{2})
                if str2double(numb{1}) + 1 == str2double(numb{2})
                    chInd1 =  numb{1};
                    if numel(chInd1) == 1
                        chInd1 = ['0' chInd1];
                    end
                    chInd2 =  numb{2};
                    if numel(chInd2) == 1
                        chInd2 = ['0' chInd2];
                    end
                    newLabel = strcat(Label(1:idx(1)-1), chInd1, '-', Label(1:idx(1)-1), chInd2);
                elseif str2double(numb{1}) == str2double(numb{2}) + 1
                    chInd1 =  numb{2};
                    if numel(chInd1) == 1
                        chInd1 = ['0' chInd1];
                    end
                    chInd2 =  numb{1};
                    if numel(chInd2) == 1
                        chInd2 = ['0'  chInd2];
                    end
                    newLabel = strcat(Label(1:idx(1)-1), chInd1, '-', Label(1:idx(1)-1), chInd2);
                end
            elseif numel(numb) == 4 && strcmp(chInd, numb{4})
                if str2double(numb{2}) + 1 == str2double(numb{4})
                    chInd1 =  numb{2};
                    if numel(chInd1) == 1
                         chInd1 = ['0' chInd1];
                    end
                    chInd2 =  numb{4};
                    if numel(chInd2) == 1
                        chInd2 = ['0' chInd2];
                    end
                    newLabel = strcat(Label(1:idx(2)-1), chInd1, '-', Label(1:idx(2)-1), chInd2);
                elseif str2double(numb{2}) == str2double(numb{4}) + 1
                    chInd1 =  numb{4};
                    if numel(chInd1) == 1
                        chInd1 = ['0' chInd1];
                    end
                    chInd2 =  numb{2};
                    if numel(chInd2) == 1
                        chInd2 = ['0' chInd2];
                    end
                    newLabel = strcat(Label(1:idx(2)-1), chInd1, '-', Label(1:idx(2)-1), chInd2);
                end 
            end            
        elseif strcmpi(montageType,'monopolar')
            if numel(chInd) == 1
                chInd = ['0'  chInd];
                newLabel = strcat(Label(1:iLastLetter), chInd);
            end
        end
        newLabels(i,1) = {newLabel};           
    end
end
        
        
