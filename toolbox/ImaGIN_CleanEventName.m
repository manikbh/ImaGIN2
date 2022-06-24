function [chLabel1, chLabel2, noteNameNew, chInd1, chInd2]  = ImaGIN_CleanEventName(Anota)
Anota = char(Anota);
xpr1  = '\w*hz_\w*';
xpr2  = '\w*stim\w*';
xpr3  = '\w*mA\w*';
xpr4  = '\w*50.0hz\w*';
xpr5  = '\w*50hz\w*';
xpr6  = '\w*50 hz\w*';

xpr4b  = '\w*55.0hz\w*';
xpr5b  = '\w*55hz\w*';
xpr6b  = '\w*55 hz\w*';

xpr7  = '\w*alarme\w*';
xpr8  = '\w*SE1Hz\w*';
xpr9  = '\w*SE 1Hz\w*';
xpr10  = 'crise';
xpr11  = '\w*OFF\w*';
xpr12  = '\w*t65535\w*';
xpr13  = '\w*Part\w*';
[~,di]=regexp(Anota,'\d*','Match');
KeepEvent = 0;
if ~isempty(di)   
    if ~isempty(regexpi(Anota,xpr7)) ||~isempty(regexpi(Anota,xpr8)) ||...
            ~isempty(regexpi(Anota,xpr9)) || ~isempty(regexpi(Anota,xpr12))||...
            ~isempty(regexpi(Anota,xpr11))|| ~isempty(regexpi(Anota,xpr13))||...
             strcmpi(Anota(1:min([length(Anota) 5])),xpr10)
    elseif ~isempty(regexpi(Anota,xpr1)) ||~isempty(regexpi(Anota,xpr4))||...
           ~isempty(regexpi(Anota,xpr5)) ||~isempty(regexpi(Anota,xpr6))||...
           ~isempty(regexpi(Anota,xpr4b))||~isempty(regexpi(Anota,xpr5b))||...
           ~isempty(regexpi(Anota,xpr6b))
        KeepEvent = 1;
    elseif ~isempty(regexpi(Anota,xpr2))
        KeepEvent = 1;
    elseif ~isempty(regexp(Anota,xpr3,'ONCE'))
        KeepEvent = 1;
    elseif ismember(lower(Anota(1)),['a':'z']) && di(1)<=4 && ~strcmp(regexprep(Anota,' ',''),'SE1Hz') && ~strcmp(regexprep(Anota,' ',''),'SE50Hz')
        KeepEvent = 1;
    end
end

if KeepEvent == 1 % Navigate all stim events
    
    % Detect stimulations and stimulation indices & save in text file
    noteName = strrep(char(Anota), ' ','_');
    noteName = regexprep(noteName,'�sec','us');
    noteName = regexprep(noteName,'�s','us');
    noteName(~ismember(double(noteName),['A':'Z' 'a':'z' '_' '.' '''' 'µ' '-' '0':'9'])) ='_';
    noteName = regexprep(noteName,'_+','_'); noteName = regexprep(noteName,'µ','u');
    noteName = strrep(noteName,'usec','us');
    noteName = strrep(noteName,'�sec','us'); % HUH notes
    noteName = regexprep(noteName,'MA','mA'); %OD
    noteName = regexprep(noteName,'Stim_Start_',''); %YUQ notes
    noteName = regexprep(noteName,'Stim_Stop_','');  %YUQ notes
    noteName = strrep(noteName,'-','_');  noteName = strrep(noteName,'__','_');
    noteName = strrep(noteName,',','');noteName = strrep(noteName,'_mA_','_');
    noteName = strrep(noteName,'sec','us');
    noteName = strrep(noteName,'_us','us');
    noteName = strrep(noteName,'AA','A'); noteName = strrep(noteName,'_MA_','_'); %some MIL notes
    keepN = ''; noteName = strrep(noteName,'stim','');  noteName = strrep(noteName,'Stim','');
    noteName = strrep(noteName,'TextNote:',''); % for BRN datasets
    
    [numZ, numZI] = regexp(noteName,'\d*','Match');
    
    try
        fundc = strfind(noteName,'_');
        lNumb = strfind(noteName,noteName(1:fundc(1)-1));
        iLastLetter = find(~ismember(noteName(1:fundc(1)-1),   '0123456789'), 1, 'last');
        keepN = noteName(1:fundc(1)-1);
        if(numel(lNumb)) == 2 && ~strcmp(keepN,'A') && ~strcmp(keepN,'H') && ~isempty(keepN)
            noteName = strrep(noteName,keepN,'CHNAME');
        elseif iLastLetter > 0 && iLastLetter > numZI(1)
            keepN = noteName(1:iLastLetter);
            noteName = strrep(noteName,keepN,'CHNAME');
        end
    end
    noteName = strrep(noteName,'.0',''); noteName = strrep(noteName,'.','');
    idScore = strfind(noteName,'_');
    if ~isempty(idScore)
        if idScore(1) < numZI(1)
            noteName(idScore(1)) = '';
        end
    end
    %% check if stim electr numbers are concatenated without space or -
    
    numbr = regexp(noteName,'\d*','Match');
    if numel(numbr) > 2 && numel(numbr) <= 5
        if str2double(numbr(2))~= str2double(numbr(1)) + 1 && str2double(numbr(1))~= str2double(numbr(2)) + 1
            if numel(numbr{1}) == 2
                if numel(numbr{2}) ~= 3 && numel(numbr{2}) ~= 2 && numel(numbr{2}) ~= 4
                    if str2double(numbr{1}(1)) + 1 == str2double(numbr{1}(2)) || str2double(numbr{1}(1)) == str2double(numbr{1}(2)) + 1
                        elecno = strcat(numbr{1}(1),'_',numbr{1}(2));
                        noteName = strrep(noteName,numbr{1},elecno);
                    end
                elseif numel(numbr{2}) == 3
                    if str2double(numbr{2}(1:2)) == str2double(numbr(1)) + 1 || str2double(numbr{2}(1:2)) + 1 == str2double(numbr(1))
                        elecno = strcat(numbr{2}(1:2),'_',numbr{2}(3));
                        noteName = regexprep(noteName,numbr{2},elecno,'once'); % We should replace the first match only
                    elseif str2double(numbr{2}(1)) == str2double(numbr(1)) + 1 || str2double(numbr{2}(1)) + 1 == str2double(numbr(1))
                        elecno = strcat(numbr{2}(1),'_',numbr{2}(2:3));
                        noteName = strrep(noteName,numbr{2},elecno);
                    end
                elseif numel(numbr{2}) == 4
                    if str2double(numbr{2}(1:2)) == str2double(numbr(1)) + 1 || str2double(numbr{2}(1:2)) + 1 == str2double(numbr(1))
                        elecno = strcat(numbr{2}(1:2),'_',numbr{2}(3:4));
                        noteName = strrep(noteName,numbr{2},elecno);
                    end
                elseif numel(numbr{2})== 2
                    if str2double(numbr{2}(1)) == str2double(numbr(1)) + 1 || str2double(numbr{2}(1)) + 1 == str2double(numbr(1))
                        elecno = strcat(numbr{2}(1),'_',numbr{2}(2));
                        noteName = strrep(noteName,numbr{2},elecno);
                    end
                end
            elseif numel(numbr{1}) == 3
                if str2double(numbr{1}(1))== 9 && str2double(numbr{1}(2:3))== 10
                    elecno = strcat(numbr{1}(1),'_',numbr{1}(2:3));
                    noteName = strrep(noteName,numbr{1},elecno);
                elseif str2double(numbr{1}(1:2))== 10 && str2double(numbr{1}(2:3))== 9
                    elecno = strcat(numbr{1}(1:2),'_',numbr{1}(3));
                    noteName = strrep(noteName,numbr{1},elecno);
                elseif str2double(numbr{1}(1)) + 1 == str2double(numbr{1}(2)) || str2double(numbr{1}(1)) == str2double(numbr{1}(2)) + 1
                    elecno = strcat(numbr{1}(1),'_',numbr{1}(2),'_',numbr{1}(3));
                    noteName = strrep(noteName,numbr{1},elecno);
                end
            elseif numel(numbr{1}) == 4
                if str2double(numbr{1}(1:2)) + 1 ==  str2double(numbr{1}(3:4)) || str2double(numbr{1}(1:2)) ==  str2double(numbr{1}(3:4)) + 1
                    elecno = strcat(numbr{1}(1:2),'_',numbr{1}(3:4));
                    noteName = strrep(noteName,numbr{1},elecno);
                elseif str2double(numbr{1}(1))== 9 && str2double(numbr{1}(2:3)) == 10
                    elecno = strcat(numbr{1}(1),'_',numbr{1}(2:3),'_',numbr{1}(4));
                    noteName = strrep(noteName,numbr{1},elecno);
                elseif str2double(numbr{1}(1:2))== 10 && str2double(numbr{1}(3)) == 9
                    elecno = strcat(numbr{1}(1:2),'_',numbr{1}(3),'_',numbr{1}(4));
                    noteName = strrep(noteName,numbr{1},elecno);
                end
            elseif numel(numbr{1}) == 5
                if str2double(numbr{1}(1:2)) + 1 ==  str2double(numbr{1}(3:4)) || str2double(numbr{1}(1:2)) ==  str2double(numbr{1}(3:4)) + 1
                    elecno = strcat(numbr{1}(1:2),'_',numbr{1}(3:4),'_',numbr{1}(5));
                    noteName = strrep(noteName,numbr{1},elecno);
                end
            elseif numel(numbr{1}) == 1
                if numel(numbr{2}) == 2
                    if str2double(numbr{2}(1)) == str2double(numbr(1)) + 1
                        elecno = strcat(numbr{2}(1),'_',numbr{2}(2));
                        noteName = strrep(noteName,numbr{2},elecno);
                    else
                        elecno = strcat(numbr{2}(1),'_',numbr{2}(2));
                        noteName = strrep(noteName,numbr{2},elecno);
                    end
                elseif numel(numbr{2}) == 3
                    if str2double(numbr{2}(1)) == str2double(numbr(1)) + 1 || str2double(numbr{2}(1)) + 1 == str2double(numbr(1))
                        elecno = strcat(numbr{2}(1),'_',numbr{2}(2:3));
                        noteName = strrep(noteName,numbr{2},elecno);
                    else
                        if str2double(numbr(1)) + 1 == str2double(numbr{2}(1:2))
                            elecno = strcat(numbr{2}(1:2),'_',numbr{2}(3));
                            noteName = strrep(noteName,numbr{2},elecno);
                        end
                    end
                end
            end
        end
    elseif numel(numbr) == 2
        if ~isempty(strfind(noteName, [numbr{2} 'us'])) || ~isempty(strfind(noteName, [numbr{2} 'mA'])) ...
                || ~isempty(strfind(noteName, [numbr{2} 'Hz']))
            if numel(numbr{1}) == 2
                elecno = strcat(numbr{1}(1),'_',numbr{1}(2));
                if str2double(numbr{1}(1)) + 1 == str2double(numbr{1}(2)) || str2double(numbr{1}(1)) == str2double(numbr{1}(2)) + 1
                    noteName = strrep(noteName,numbr{1},elecno);
                end
            elseif ~isempty(strfind(numbr{1},'10'))
                noteName = strrep(noteName,'10','_10_');
            elseif numel(numbr{1}) == 4
                elecno = strcat(numbr{1}(1:2),'_',numbr{1}(3));
                noteName = strrep(noteName,numbr{1},elecno);
            elseif ~isempty(strfind(noteName, [numbr{1} 'mA'])) || ~isempty(strfind(noteName, [numbr{1} 'Hz']))
                if isempty(strfind(numbr{1},'.'))
                    if numel(numbr{1}) == 3
                        elecno = strcat(numbr{1}(1),'_',numbr{1}(2),'_',numbr{1}(3));
                        noteName = strrep(noteName,numbr{1},elecno);
                    elseif numel(numbr{1}) == 5
                        elecno = strcat(numbr{1}(1:2),'_',numbr{1}(3:4),'_',numbr{1}(5));
                        noteName = strrep(noteName,numbr{1},elecno);
                    end
                end
            end
        end
    end
    
    numb = regexp(noteName,'\d*','Match');
    
    if numel(numb) >= 2
        idx1 = strfind(noteName,numb{1});
        idx2 = strfind(noteName,numb{2});
        cnbre1= numel(numb{1});
        cnbre2= numel(numb{2});
        if str2double(numb(1)) == str2double(numb(2)) - 1
            sfix = noteName(1:(idx1(1)+cnbre1)-1);
            noteName = strcat(upper(sfix),noteName(idx2(1):end));
        elseif str2double(numb(1)) == str2double(numb(2)) + 1
            sfix = noteName(1:idx1(1)-1);
            noteName = strcat(upper(sfix),num2str(numb{2}), num2str(numb{1}),noteName(idx2(1)+cnbre2:end));
        else
            if str2double(numb(1)) > str2double(numb(2))
                sfix = noteName(1:idx1(1)-1);
                noteName = strcat(upper(sfix),num2str(numb{2}), num2str(numb{1}),noteName(idx2(1)+cnbre2:end));
            elseif str2double(numb(1)) < str2double(numb(2))
                sfix = noteName(1:idx1(1)-1);
                noteName = strcat(upper(sfix),num2str(numb{1}), num2str(numb{2}),noteName(idx2(1)+cnbre2:end));
            end
        end
        
        xpr1  = '\w*Hz_\w*'; xpr2 = '\w*us_\w*';
        xpri1 = regexpi(noteName,xpr1); xpri2 = regexpi(noteName,xpr2);
    end
    noteName = strrep(noteName,'-','');  noteName = strrep(noteName,'_mA','mA');
    
    %% build .mat/.dat name
    
    try
        [ds,di] = regexp(noteName,'\d*','Match');
        xsub0 = noteName(1:(di(1)+numel(ds{1})-1));
        rxp1  = '[-+]?(\d*[.])?\d+mA'; rxp2  = '[-+]?(\d*[.])?\d+Hz';
        rxp3  = '[-+]?(\d*[.])?\d+us'; rxp4  = '[-+]?(\d*[.])?\d+s';
        xsub1 = regexp(noteName,rxp1,'match');
        if isempty(xsub1), xsub1 = '0mA';
            xsub1 = cellstr(xsub1);
        end
        xsub2 = regexp(noteName,rxp2,'match');
        if isempty(xsub2), xsub2 = '0Hz';
            xsub2 = cellstr(xsub2);
        end
        xsub3 = regexp(noteName,rxp3,'match');
        xsub4 = regexp(noteName,rxp4,'match');
        if ~isempty(xsub3)&& isempty(xsub4)
            xsub3 = cellstr(xsub3);
        elseif isempty(xsub3)&& ~isempty(xsub4)
            xsub4 = char(xsub4);
            if numel(xsub4) <= 3
                buff = num2str(str2double(xsub4));
                if numel(buff) <= 4
                    xsub4 = buff;
                end
            end
            xsub3 = cellstr(strcat(xsub4(1:end-1),'us'));
        else
            xsub3 = '0us';
        end
        
        FullN = strcat(upper(xsub0),'_',xsub1,'_',xsub2,'_',xsub3);
        FullN = char(unique(FullN));
        FullN = FullN(1,:);
        if isempty(FullN)
            FullN = noteName;
        end
    catch
        FullN = noteName;
    end
    
    noteName = strrep(FullN,'.',',');
    
    %%
    numb2 = regexp(noteName,'\d*','Match');
    if numel(numb2) >= 1
        idxn1 = strfind(noteName,numb2{1});          %OD
        subn1 = strrep(noteName(1:idxn1),'_','');
        noteName = char(strcat(subn1,noteName(idxn1+1:end)));
    end
    
    noteName = strrep(noteName,'CHNAME',upper(keepN));
    
    ptrn = ',';
    if strncmp(noteName,ptrn,1)
        noteName = char(noteName(2:end));
    end

    %%
    [~,tmpdi] = regexp(noteName,'\d*','Match');
    noteNameNew = noteName;
    noteNameNew(1:tmpdi(1)-1) = upper(noteNameNew(1:tmpdi(1)-1));
    idxScore = strfind(noteNameNew,'_');
    if isempty(idxScore)
        chInd1 = '';
        chInd2 = '';
        chLabel1 = Anota;
        chLabel2 = Anota;
        noteNameNew = Anota;
        return;
    end
    Label = noteNameNew(1:idxScore(1)-1);
    iLastLetter = find(~ismember(Label, '0123456789'), 1, 'last');
    if isempty(iLastLetter) || (iLastLetter == length(Label))
        chInd1 = '';
        chInd2 = '';
        chLabel1 = Anota;
        chLabel2 = Anota;
        noteNameNew = Anota;
        return;        
    end
    chLabel = Label(1:iLastLetter);
    chInd = Label(iLastLetter+1:end);
    if numel(chInd)==2        
        chInd1 = chInd(1);
        chInd2 = chInd(2);
        chLabel1 = strcat(chLabel, chInd1);
        chLabel2 = strcat(chLabel, chInd2);       
        noteNameNew = strcat(chLabel1, '-',chLabel2, noteNameNew(idxScore(1):end));
    elseif numel(chInd)==3
        chInd1 = chInd(1);
        chLabel1 = strcat(chLabel, chInd1);        
        chInd2 = chInd(2:3);
        chLabel2 = strcat(chLabel,  chInd2);
        noteNameNew = strcat(chLabel1,'-',chLabel2, noteNameNew(idxScore(1):end));
    elseif numel(chInd)==4
        chInd1 = chInd(1:2);
        chInd2 = chInd(3:4);
        chLabel1 = strcat(chLabel,chInd1);
        chLabel2 = strcat(chLabel, chInd2);
        noteNameNew = strcat(chLabel1, '-',chLabel2, noteNameNew(idxScore(1):end));
    elseif numel(chInd)==6
        chInd1 = chInd(1:3);
        chInd2 = chInd(4:6);
        chLabel1 = strcat(chLabel, chInd1);
        chLabel2 = strcat(chLabel, chInd2);
        noteNameNew = strcat(chLabel1, '-',chLabel2, noteNameNew(idxScore(1):end));
    else
        chInd1 = '';
        chInd2 = '';
        chLabel1 = Anota;
        chLabel2 = Anota;
        noteNameNew = Anota;
    end
else
   chInd1 = '';
   chInd2 = '';
   chLabel1 = Anota;
   chLabel2 = Anota;
   noteNameNew = Anota;
end