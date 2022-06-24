function [Stimulation,StimulationIndex,StimulationFreqU]=ImaGIN_StimDetect(S)

try
    Filename=S.Fname;
catch
    Filename = spm_select(1, '\.mat$', 'Select data file');
end

D=spm_eeg_load(Filename);

try
    SelChan=S.Channels;
    if isempty(SelChan)
        SelChan=setdiff(setdiff(1:nchannels(D),indchantype(D,'ECG')),badchannels(D));
    end
catch
    SelChan = spm_input('Select channel', 1, 'i');
end

try
    Start=S.StimStart;
catch
    Start=[];
end

try
    End=S.StimEnd;
catch
    End=[];
end

try
    Stim=S.StimFreq;
catch
    Stim=spm_input('Stimulation frequency [Hz]', '+1', 'r');
end

try
    minStim=S.minStim;
catch
    minStim=spm_input('Minimal number of stimulations', '+1', 'r');
end

try
    FindBadChannels=S.FindBadChannels;
catch
    FindBadChannels=1;
end

try
    StimContinuous=S.StimContinuous;
catch
    StimContinuous=spm_input('Continuous stimulation ', '+1', 'Yes|No');
    switch StimContinuous
        case 'Yes'
            StimContinuous=1;
        case 'No'
            StimContinuous=0;
    end
end

try
    FileOut=S.FileOut;
catch
    FileOut = [Filename(1:end-3) 'txt'];
end

try
    StimDetectVersion=S.StimDetectVersion;
catch
    StimDetectVersion=4;
end


% D=timeonset(D,0);
% save(D)
Time=time(D);


StimInit=Stim;
Stim=D.fsample/StimInit;
if isempty(Start)
    Start=1;
else
    Start=min(indsample(D,Start));
    if isnan(Start)
        Start=1;
    end
end
if isempty(End)
    End=nsamples(D);
else
    End=min(indsample(D,End));
    if isnan(End)
        End=nsamples(D);
    end
end

% Data=sum(D(SelChan,Start:End),1);
Data=D(SelChan,Start:End);

%find bad channels which saturate
if FindBadChannels
    L=zeros(1,size(Data,1));
    for i1=1:size(Data,1)
        tmp=abs(Data(i1,1:2:end));
        L(i1)=length(find(tmp==max(tmp)))/length(tmp);
    end
%     GoodChannels=find(L<200/size(Data,2));
    GoodChannels=find(L<0.05);      %Assume that bad channels saturate 5% of time
    Data=ImaGIN_normalisation(Data(GoodChannels,:),2);
end

%Remove line noise at 50 Hz
for i1=1:size(Data,1)
    Data(i1,:)=ImaGIN_notch(Data(i1,:),D.fsample,50,250);
end

%Remove low frequency drift (below 3 Hz)
for i1=1:size(Data,1)
    Data(i1,:)=ImaGIN_bandpass(Data(i1,:),D.fsample,1,0.95*D.fsample/2);
end


%Simple version based on averaging between electrodes
stimulation = [];
if 1==1
    
    % d=[0 abs(diff(Data,2)) 0];
    d=[zeros(size(Data,1),1) abs(diff(Data,2,2)) zeros(size(Data,1),1)];
%     d=ImaGIN_normalisation(d,2);
%     [tmp,Order]=sort(max(d'),'descend');
%     d=mean(d(Order(1:ceil(length(Order)/2)),:));


    d=ImaGIN_normalisation(d,2);
    d1=max(d,[],1);
    %to avoid outliers we take the 95 percentile instead of max as the line
    %above
    d=sort(d,1);
    npercentile=ceil(0.05*size(d,1));
    d1=d(npercentile,:);
    
    
    d1=ImaGIN_normalisation(d1,2);
    d2=mean(d,1);
    d2=ImaGIN_normalisation(d2,2);  
    if isempty(d1)  % fix: matlab2015 doesn't know how to make addition NaN matrix and empty matrix
        d1 = NaN(size(d2));
    end
    
    d=d1+d2;
    
    
    d=ImaGIN_normalisation(d,2);
    
    % [tmp1,tmp2]=max(d);
    % Index=find(d>tmp1/4);
    Th=2;
    ok=1;
    while ok
        EstimatedStim=(length(d)-ceil(Stim))./ceil(Stim);
        Index=find(d>Th);
        Index=Index(find(Index>ceil(Stim/2)+1&Index<length(d)-ceil(Stim/2)-1));
        if length(Index)/5>EstimatedStim
            Th=Th+1;
        else
            ok=0;
        end
        if Th > 1000 % VT: control infinite loop 
            ok = 0;
        end
    end
    if ~isempty(Index)
        IndexClust=zeros(size(Index));
        IndexClust(1)=1;
        for i1=2:length(Index)
            tmp1=find(abs(Index(i1)-Index)<ceil(Stim/2));
            if max(IndexClust(tmp1))==0
                IndexClust(tmp1)=max(IndexClust(1:i1-1))+1;
            else
                IndexClust(tmp1)=max(IndexClust(tmp1));
            end
            %         IndexClust(tmp1(find(IndexClust(tmp1)==0)))=min(IndexClust(1:i1-1))+1;
        end
        IndexNew=zeros(1,max(IndexClust));
        for i1=1:length(IndexNew)
            tmp=find(IndexClust==i1);
            IndexNew(i1)=Index(tmp(floor(median(1:length(tmp)))));
        end
        Index=IndexNew;
        
        Template=0;
        for i1=1:length(Index)
            Template=Template+d(Index(i1)+[-ceil(Stim/2):ceil(Stim/2)]);%./sum(abs(d(Index(i1)+[-ceil(Stim/2):ceil(Stim/2)])));
        end
        
        % Template=d(tmp2+[-ceil(Stim/2):ceil(Stim/2)]);
        
        
        cc=zeros(size(d));
        % for i1=ceil(Stim/2)+1:length(d)-ceil(Stim/2)-1
        for i1=Index
            tmp=corrcoef(Template,d(i1+[-ceil(Stim/2):ceil(Stim/2)]));
            cc(i1)=tmp(2);
            tmp1=find(i1-Index<ceil(Stim/3)&i1-Index>=0);
            if length(tmp1)>1
                cc(Index(tmp1(find(cc(Index(tmp1))<max(cc(Index(tmp1)))))))=0;
            end
        end
        %Threshold cc
        Threshold=0.4;
        Threshold=min([0.4 median(cc(find(cc~=0)))]);
        tmp=find(cc>Threshold);
        stimulation=tmp;
        
        
        %Second pass
        Index=stimulation;
        Template=0;
        for i1=1:length(Index)
            if Index(i1)-ceil(Stim/2) >= 1 && Index(i1)+ceil(Stim/2) <= length(d)
%                 Template=Template+d(Index(i1)+[-ceil(Stim/2):ceil(Stim/2)]);
                Template=Template+d(Index(i1)+[-ceil(Stim/2):ceil(Stim/2)])./sum(abs(d(Index(i1)+[-ceil(Stim/2):ceil(Stim/2)])));
            end
        end
        Index=find(d>2);
        Index=Index(find(Index>ceil(Stim/2)+1&Index<length(d)-ceil(Stim/2)-1));
        cc=zeros(size(d));
        % for i1=ceil(Stim/2)+1:length(d)-ceil(Stim/2)-1
        for i1=Index
            tmp=corrcoef(Template,d(i1+[-ceil(Stim/2):ceil(Stim/2)]));
            cc(i1)=tmp(2);
            tmp1=find(i1-Index<ceil(Stim/3)&i1-Index>=0);
            if length(tmp1)>1
                cc(Index(tmp1(find(cc(Index(tmp1))<max(cc(Index(tmp1)))))))=0;
            end
        end
        %Threshold cc
        Threshold=0.5;
        tmp=find(cc>Threshold);
        stimulation=tmp;

    end 
    
    % remove outliers close to other stims
%     stimulation = removeOutliers( stimulation, Stim ) ;
%     if StimContinuous
%             if length(stimulation)>1
%                 Stim=median(diff(stimulation)); %estimate actual stimulation frequency
%             end
            stimulation = removeOutliers( stimulation, Stim , minStim) ;
            stimulation = stimulation(find(stimulation>round( 0.95 * Stim/2 )&stimulation+round( 0.95 * Stim/2 )<=length(d)))
%     end
       
    if ~isempty(stimulation)
        
        % alignment of templates
        % 95% of segmentHWidth is enough for alignment (100% may fail when for instance, segmentHWidth = 500, first stimulation is at sample 502 and alignment is found to be -2)
        % and same pb also happens for offset correction
        % segmentHWidth       = ceil(Stim/2) ;
        segmentHWidth       = round( 0.95 * Stim/2 ) ;
        templateHWidth      = round( 0.01 * fsample(D) );
        corrHWidth          = round( 0.02 * fsample(D) ) ;
        do_plot             = false ;
        iterMax             = 10 ; 
        if segmentHWidth > templateHWidth % sometimes it happens to be not the case: 0007BUC/ElectrodeFile#2/S.EvtName = P10-P11_2mA_1011Hz_3000us
            
            [ alignments, alignedStimulations ]   = ImaGIN_AlignSegmentsIter( d, stimulation, segmentHWidth, templateHWidth, corrHWidth, do_plot, iterMax ) ;
            
            if 0
                [ pathstr, name ]   = fileparts(FileOut) ;
                fname = fullfile( pathstr, 'templateAlignment', [ name '__templateAlignment.mat' ] ) ;
                if ~exist( fileparts(fname), 'dir' ), mkdir( fileparts(fname) ) ; end
                fprintf( 1, [ 'Save ' fname '\n' ] ) ;
                save( fname, 'stimulation', 'alignments', 'alignedStimulations' ) ;
            end
            
            % apply alignment only if convergence and max(abs(alignments)) < 10ms
            for i=1:iterMax
                % convergence
                if isempty(find(alignments(:,i)~=0,1))
                    % sanity check
                    if ~isequal( alignedStimulations, stimulation+sum(alignments(:,1:i),2)' ), error( 'Inconsistent alignment' ) ; end
                    % apply alignment if < 10ms
                    if max(abs(alignedStimulations-stimulation))/fsample(D) < 10e-3
                        stimulation = alignedStimulations ;
                    end
                    break
                end
            end
        end
        
        
        % try to correct onset inpired from Lena's method
        tmps=sign(mean(Data(:,stimulation),2));
        Data2=Data;
        for i1=1:size(Data2,1)
            Data2(i1,:)=tmps(i1)*Data2(i1,:);
        end
        
        TemplateSamples=-segmentHWidth:segmentHWidth ;
        TimeTemplate=TemplateSamples./fsample(D);
        avg_to_plot = 0 ;
        for i1=1:length(stimulation)
            avg_to_plot = avg_to_plot + Data(:,stimulation(i1)+TemplateSamples) ;
        end
        
        %Same offset for all pulses
%         StimDetectVersion = 1 ; % based on template of data + threshold from 10 to 5
%         StimDetectVersion = 2 ; % based on template of acceleration d + threshold (perc of the max)
%         StimDetectVersion = 3 ; % mixed method: accel used to detect the start of the stim and then data to detect the max (threshold from 10 to 4)
%         StimDetectVersion = 4 ; % similar to 3 but offset is computed from data based on a percentage of the max
        if StimDetectVersion == 1 
            
            TemplateData=0;
            for i1=1:length(stimulation)
                TemplateData=TemplateData+mean(abs(Data2(:,stimulation(i1)+TemplateSamples)),1);
            end
            TemplateNorm=ImaGIN_normalisation(TemplateData,2,find(TimeTemplate<-0.02&TimeTemplate>-0.1));
            template_to_plot = TemplateNorm ;
            
            for th=10:-1:5
                tmpoffset=find( abs(TemplateNorm)>th & abs(TimeTemplate)<0.02);
                if ~isempty(tmpoffset)
                    break
                end
            end
            if isempty(tmpoffset)
                offset=nan;
            else
                i1=min(tmpoffset);
                offset=ceil(Stim/2)+1-i1; 
            end
            
            
        elseif StimDetectVersion == 2
      
            % use mean template from d instead of data
            template = 0 ;
            for s=1:length(stimulation)
                template = template + d( stimulation(s) + (-segmentHWidth:segmentHWidth) ) ;
            end
            template_to_plot = template ;
            
            % use a percentage of the max as a threshold
            % sometimes there are many peaks (biphasic, high fsample) and the first is not the biggest one
            % so we prefer a low threshold on the max
            th              = 0.3 ;
            offsetMaxTime   = 0.03 ;
            offsetMaxSample = round( offsetMaxTime * fsample(D) ) ;
            segmentCenter   = segmentHWidth + 1 ;
            [ templateMax, iMax ] = max( abs(template(segmentCenter+(-offsetMaxSample:offsetMaxSample))) ) ;
            
            tmpoffset       = find( template>th*templateMax & abs(TimeTemplate)<offsetMaxTime, 1 );
            
            if isempty(tmpoffset)
                offset=nan;
            else
                i1=min(tmpoffset);
                offset=segmentCenter-i1;
            end
            
            
        elseif StimDetectVersion == 3
            
            % offset is computed from the data
            TemplateData=0;
            for i1=1:length(stimulation)
                TemplateData=TemplateData+mean(abs(Data2(:,stimulation(i1)+TemplateSamples)),1);
            end
            TemplateNorm=ImaGIN_normalisation(TemplateData,2,find(TimeTemplate<-0.02&TimeTemplate>-0.1));
            
            % try to detect stimulation start to refine detection
            % use mean template from d instead of data
            template = 0 ;
            for s=1:length(stimulation)
                template = template + d( stimulation(s) + TemplateSamples ) ;
            end
            % use a percentage of the max as a threshold
            % sometimes there are many peaks (biphasic, high fsample) and the first is not the biggest one
            % so we prefer a low threshold on the max
            th                      = 30/100 ;
            offsetMaxInterval       = [ -0.04 0.02 ];
            offsetMaxIntervalSample = round( offsetMaxInterval * fsample(D) ) ;
            segmentCenter           = segmentHWidth + 1 ;
            [ templateMax, ~ ]      = max( abs(template(segmentCenter+(offsetMaxIntervalSample(1):offsetMaxIntervalSample(2)))) ) ;
            
            startoffset             = find( template>th*templateMax & TimeTemplate>offsetMaxInterval(1) & TimeTemplate<offsetMaxInterval(2), 1 );
            
            if isempty(startoffset)
                
                for th=10:-1:4
                    tmpoffset=find( abs(TemplateNorm)>th & abs(TimeTemplate)<0.02);
                    if ~isempty(tmpoffset)
                        break
                    end
                end
                
            else
                
                startoffset=min(startoffset);
                for th=10:-1:4
                    tmpoffset=find( abs(TemplateNorm)>th & TimeTemplate>=TimeTemplate(startoffset+1) & TimeTemplate<=TimeTemplate(startoffset)+8e-3 );
                    if ~isempty(tmpoffset)
                        break
                    end
                end
                if isempty(tmpoffset)
                    % 4ms after startoffset (one sample at 256Hz)
                    tmpoffset = startoffset + round( fsample(D)*4e-3 ) ;
                end
            end
            
            if isempty(tmpoffset)
                offset=nan;
            else
                i1=min(tmpoffset);
%                 offset=ceil(Stim/2)+1-i1;
                offset=segmentCenter-i1;
            end
            
            
        elseif StimDetectVersion == 4
            
            % offset is computed from the data
            TemplateData=0;
            for i1=1:length(stimulation)
                TemplateData=TemplateData+mean(abs(Data2(:,stimulation(i1)+TemplateSamples)),1);
            end
            TemplateNorm=ImaGIN_normalisation(TemplateData,2,find(TimeTemplate<-0.02&TimeTemplate>-0.1));
            
            % try to detect stimulation start to refine detection
            % use mean template from d instead of data
            template = 0 ;
            for s=1:length(stimulation)
                template = template + d( stimulation(s) + TemplateSamples ) ;
            end
            % use a percentage of the max as a threshold
            % sometimes there are many peaks (biphasic, high fsample) and the first is not the biggest one
            % so we prefer a low threshold on the max
            th                      = 30/100 ;
            offsetMaxInterval       = [ -0.04 0.02 ];
            offsetMaxIntervalSample = round( offsetMaxInterval * fsample(D) ) ;
            % [ -0.04 0.02 ] is adapted only for 1Hz stimulation, not for 50Hz
            % so we restrict the interval to avoid errors occuring in case offsetMaxInterval is larger than TemplateSamples
            offsetMaxIntervalSample = [ max(offsetMaxIntervalSample(1),TemplateSamples(1)) min(offsetMaxIntervalSample(2),TemplateSamples(end)) ];
            segmentCenter           = segmentHWidth + 1 ;
            [ templateMax, ~ ]      = max( abs(template(segmentCenter+(offsetMaxIntervalSample(1):offsetMaxIntervalSample(2)))) ) ;
            
            startoffset             = find( template>th*templateMax & TimeTemplate>offsetMaxInterval(1) & TimeTemplate<offsetMaxInterval(2), 1 );
            
            if isempty(startoffset)
                
                for th=10:-1:4
                    tmpoffset=find( abs(TemplateNorm)>th & abs(TimeTemplate)<0.02);
                    if ~isempty(tmpoffset)
                        break
                    end
                end
                
            else
                
                startoffset = min(startoffset);
                % 6ms allow to move 2 samples away at 256 Hz
                % when the Note is bizarre, the Stim frequency is not accurate and TemplateSamples is not well defined (only a few samples) 
                % so we need to guarantee that startoffset + 6ms do not fall outside TemplateSamples range
                tmpEndSamp  = min(startoffset+round(fsample(D)*6e-3),length(TemplateNorm)) ;
                tmpMax   	= max( abs(TemplateNorm(startoffset:tmpEndSamp)) ) ;
                thd         = 50/100 ;
                tmpoffset   = find( abs(TemplateNorm)>thd*tmpMax & TimeTemplate>=TimeTemplate(startoffset) & TimeTemplate<=TimeTemplate(tmpEndSamp) );
                % we take the sample at the first third of the interval (usefull for bipolar stim at 4096Hz)
                tmpoffset   = tmpoffset( floor((length(tmpoffset)-1)/3)+1 ) ;
            end
            
            if isempty(tmpoffset)
                offset=nan;
            else
                i1=min(tmpoffset);
%                 offset=ceil(Stim/2)+1-i1;
                offset=segmentCenter-i1;
            end            
            
        elseif 1==1
            %Same offset for all pulses
            TemplateData=0;
            for i1=1:length(stimulation)
%                 TemplateData=TemplateData+mean(Data2(:,stimulation(i1)+[-ceil(Stim/2):ceil(Stim/2)]),1);
                TemplateData=TemplateData+mean(abs(hilbert(Data2(:,stimulation(i1)+[-ceil(Stim/2):ceil(Stim/2)]))),1);
%                 TemplateData=TemplateData+mean(abs(Data2(:,stimulation(i1)+[-ceil(Stim/2):ceil(Stim/2)])),1);
%                 TemplateData=TemplateData+mean(ImaGIN_normalisation(abs(hilbert(Data2(:,stimulation(i1)+[-ceil(Stim/2):ceil(Stim/2)]))),2,find(TimeTemplate<-0.02&TimeTemplate>-0.1)),1);
            end
            TimeTemplate=[-ceil(Stim/2):ceil(Stim/2)]./fsample(D);
            TemplateNorm=ImaGIN_normalisation(TemplateData,2,find(TimeTemplate<-0.01&TimeTemplate>-0.1));
%             tmpoffset=find(abs(hilbert(TemplateNorm))>10&TimeTemplate>-0.02);
            tmpoffset=find(abs(TemplateNorm)>10&TimeTemplate>-0.02);
            if isempty(tmpoffset)
                offset=0;
            else
                i1=min(tmpoffset);
                offset=ceil(Stim/2)+1-i1;
            end
        else
            %different offset for all pulses
            TimeTemplate=[-ceil(Stim/2):ceil(Stim/2)]./fsample(D);
            offset=stimulation;
            for i1=1:length(stimulation)
                TemplateData=mean(Data2(:,stimulation(i1)+[-ceil(Stim/2):ceil(Stim/2)]),1);
                TemplateNorm=ImaGIN_normalisation(TemplateData,2,find(TimeTemplate<-0.01&TimeTemplate>-0.1));
                tmpoffset=find(abs(hilbert(TemplateNorm))>10&TimeTemplate>-0.015);
                if isempty(tmpoffset)
                    offset(i1)=0;
                else
                    i2=min(tmpoffset);
                    offset(i1)=ceil(Stim/2)+1-i2;
                end
%                 if abs(offset(i1))>5
%                     offset(i1)=0;
%                 end                    
            end
        end
        
        
        % save and plot for comparison of methods
        plot_and_save = false ;
        if plot_and_save
            if StimDetectVersion == 1 || StimDetectVersion == 2 || StimDetectVersion == 3 || StimDetectVersion == 4
                
                % save th and offset in txt file
                [ pathstr, name ] = fileparts(FileOut) ;
                fname = fullfile( pathstr, 'offset', [ name '_offset.txt' ] ) ;
                if ~exist(fileparts(fname),'dir'), mkdir(fileparts(fname)); end
                fprintf( 1, [ 'Save ' fname '\n' ] ) ;
                fd = fopen( fname, 'w' ) ;
                fprintf( fd, [ num2str(offset) ' ' num2str(offset/fsample(D)) '\n' ] ) ;
                fprintf( fd, [ num2str(th) ] ) ;
                fclose(fd);
                
                % plot offset
                H = figure( 'Position', [189 118 560 700] ) ; %[30 183 560 700] ) ;
                
                subplot(3,1,1); hold on ;
                plot( TimeTemplate, avg_to_plot ) ;
                line( zeros(2,1), get(gca,'YLim')','Color','k','linewidth',2);
                if ~isnan(offset)
                    line([TimeTemplate(i1);TimeTemplate(i1)],get(gca,'YLim')','color','r','linewidth',2);
                end
                xlim( [ -.05 .05 ] ) ;
                title( strrep( [ name ' - V' num2str(StimDetectVersion) ], '_', ' ' ) )
                
                subplot(3,1,2); hold on ;
                plot( TimeTemplate, template,'Color','b','linewidth',3 ) ;
                line( zeros(2,1), get(gca,'YLim')','Color','k','linewidth',2);
                if ~isempty(startoffset)
                    line([TimeTemplate(startoffset);TimeTemplate(startoffset)],get(gca,'YLim')','color','b','linewidth',2);
                end
                xlim( [ -.05 .05 ] ) ;
                title( 'template accel' )
                
                subplot(3,1,3); hold on ;
                plot( TimeTemplate, TemplateNorm,'Color','g','linewidth',3 ) ;
                line( zeros(2,1), get(gca,'YLim')','Color','k','linewidth',2);
                if ~isnan(offset)
                    line([TimeTemplate(i1);TimeTemplate(i1)],get(gca,'YLim')','color','g','linewidth',2);
                end
                xlim( [ -.05 .05 ] ) ;
                title( 'template data' )
                
                fname = fullfile( pathstr, 'offset', [ name '_offset.png' ] ) ;
                fprintf( 1, [ 'Save ' fname '...' ] ) ;
                saveas(H,fname);
                fprintf( 1, [ 'done' '\n' ] ) ;

                close(H);
                
            end
        end
        
        
        if ~isnan(offset)
            stimulation=stimulation-offset;
        end
    end
    
else
    
    %more complex version where start and end stim are automatically found
    
    clear stimulation StartStim EndStim template
    n=0;
    for i0=1:size(Data,1)
        % d=[0 abs(diff(Data,2)) 0];
        d=[0 abs(diff(Data(i0,:),2,2)) 0];
        d=ImaGIN_normalisation(d,2);
        
        % [tmp1,tmp2]=max(d);
        % Index=find(d>tmp1/4);
        Index=find(d>3);
        
        Index=Index(find(Index>ceil(Stim/2)+1&Index<length(d)-ceil(Stim/2)-1));
        Template=0;
        for i1=1:length(Index)
            Template=Template+d(Index(i1)+[-ceil(Stim/2):ceil(Stim/2)]);
        end
        
        % Template=d(tmp2+[-ceil(Stim/2):ceil(Stim/2)]);
        
        
        cc=zeros(size(d));
        % for i1=ceil(Stim/2)+1:length(d)-ceil(Stim/2)-1
        for i1=Index
            tmp=corrcoef(Template,d(i1+[-ceil(Stim/2):ceil(Stim/2)]));
            cc(i1)=tmp(2);
            tmp1=find(i1-Index<ceil(Stim/3)&i1-Index>=0);
            if length(tmp1)>1
                cc(Index(tmp1(find(cc(Index(tmp1))<max(cc(Index(tmp1)))))))=0;
            end
        end
        
        %Threshold cc
        Threshold=0.25;
        tmp=find(cc>Threshold);
        stimulation{i0}=tmp;
        if ~isempty(tmp)
            StartStim(i0)=min(tmp);
            EndStim(i0)=max(tmp);
            n=n+1;
            template(n,:)=Template;
        end
        
    end
    
    Template=mean(template);
    CC=zeros(1,size(template,1));
    for i0=1:length(CC)
        cc=corrcoef(template(i0,:),Template);
        CC(i0)=cc(2);
    end
    Index=find(StartStim>0);
    Index=Index(find(CC>0.6));
    
        
% % %     StartStim=median(StartStim(find(EndStim>0)))-fsample(D);
% % %     EndStim=median(EndStim(find(EndStim>0)))+fsample(D);
% %     StartStim=min(StartStim(find(EndStim>0)))-fsample(D);
% %     EndStim=max(EndStim(find(EndStim>0)))+fsample(D);
%     
%     
%     Duration=(EndStim-StartStim)./fsample(D);
%     Index=find(StartStim>0);
% %     idx=kmeans(Duration(Index),2);
%     idx=kmeans([StartStim(Index);EndStim(Index);Duration(Index)]',2);
%     if mean(Duration(Index(find(idx==1))))>mean(Duration(Index(find(idx==2))))
%         Index=Index(find(idx==1));
%     else
%         Index=Index(find(idx==2));
%     end
%     
    
    
    StartStim=min(StartStim(Index))-fsample(D);
    EndStim=max(EndStim(Index))+fsample(D);
    
    
    d=[zeros(size(Data(Index,:),1),1) abs(diff(Data(Index,:),2,2)) zeros(size(Data(Index,:),1),1)];
    % d=ImaGIN_normalisation(d,2);
    % [tmp,Order]=sort(max(d'),'descend');
    % d=mean(d(Order(1:ceil(length(Order)/2)),:));
    d=mean(d);
    d=ImaGIN_normalisation(d,2);
    
    % [tmp1,tmp2]=max(d);
    % Index=find(d>tmp1/4);
    Index=find(d>3);
    
    Index=Index(find(Index>=StartStim&Index<=EndStim));
    
    Index=Index(find(Index>ceil(Stim/2)+1&Index<length(d)-ceil(Stim/2)-1));
    Template=0;
    for i1=1:length(Index)
        Template=Template+d(Index(i1)+[-ceil(Stim/2):ceil(Stim/2)]);
    end
    
    % Template=d(tmp2+[-ceil(Stim/2):ceil(Stim/2)]);
    
    
    cc=zeros(size(d));
    % for i1=ceil(Stim/2)+1:length(d)-ceil(Stim/2)-1
    for i1=Index
        tmp=corrcoef(Template,d(i1+[-ceil(Stim/2):ceil(Stim/2)]));
        cc(i1)=tmp(2);
        tmp1=find(i1-Index<ceil(Stim/3)&i1-Index>=0);
        if length(tmp1)>1
            cc(Index(tmp1(find(cc(Index(tmp1))<max(cc(Index(tmp1)))))))=0;
        end
    end
    
    %Threshold cc
    Threshold=0.4;
    tmp=find(cc>Threshold);
    stimulation=tmp;
        
end


%remove outliers close to other stims or isolated clusters
stimulation = removeOutliers( stimulation, Stim , minStim) ;


if isempty(stimulation)
    StimulationFreqU=[];
    StimulationIndex = [];
    Stimulation = [];
else
    Index=find(d(stimulation(1):stimulation(end))>2);
    Index=Index(find(Index>ceil(Stim/2)+1&Index<length(d)-ceil(Stim/2)-1));
    if length(Index>1) > 2 % Index should be great than 2 not 0
        tmp=diff(Index);
        tmp=tmp(find(tmp>max([2 D.fsample/100])));
        if ~isempty(tmp)
            tmp=sort(tmp);
            tmp=tmp(ceil(length(tmp)/3):end);
            StimulationFreqU=D.fsample/median(tmp);   %uncorrected stim frequency
        else
            StimulationFreqU = 0;
        end
    else
        StimulationFreqU=0;   %uncorrected stim frequency
    end
    
    
    
    %fill the gaps
    if StimContinuous
        tmp=stimulation;
        stimulation=tmp(1);
        while stimulation(end)<tmp(end)
            [tmp1 tmp2]=min(abs(tmp-(stimulation(end)+Stim)));
            if tmp1>Stim/2
                stimulation=[stimulation stimulation(end)+Stim];
            else
                stimulation=[stimulation tmp(tmp2)];
            end
        end
    end
    
    StimulationIndex=stimulation+Start-1;
    Stimulation=Time(round(StimulationIndex));
end

% 
% tmp2=diff(tmp);
% tmp3=[1 find(tmp2~=1)+1 length(tmp)];
% Stimulation=zeros(length(tmp3)-1,1);
% for i1=1:length(tmp3)-1
%     [t1,t2]=max(cc(tmp(tmp3(i1):tmp3(i1+1))));
%     Stimulation(i1)=tmp(t2+tmp3(i1)-1);
% end
% Stimulation=Time(Start-1+Stimulation);
% 
% ok=1;
% while ok
%     StimulationOK=ones(size(Stimulation));
%     for i1=2:length(Stimulation)
%         if Stimulation(i1)-Stimulation(i1-1)<1/(10*StimInit)
%             StimulationOK(i1)=0;
%         end
%     end
%     Index=find(StimulationOK==0);
%     if isempty(Index)
%         ok=0;
%     else
%         Stimulation=Stimulation(find(StimulationOK));
%     end
% end


% P = spm_str_manip(Filename, 'H');
% F = spm_str_manip(Filename, 'tr');

fid=fopen(FileOut,'w');
fprintf(fid,'%f\n',Stimulation);
fclose(fid);

% Filename=fullfile(DirOut,[F '_StimulationIndex.txt']);
% fid=fopen(Filename,'w');
% fprintf(fid,'%f\n',StimulationIndex);
% fclose(fid);

end


%remove outliers close to other stims or isolated
function stimulation = removeOutliers( stimulation, Stim , minStim)

if length(stimulation)>1
    remove=[];
    for i1=1:length(stimulation)
        difftmp=abs(stimulation-stimulation(i1));
        if isempty(find(difftmp>0.6*Stim & difftmp<1.4*Stim,1))
            remove=[remove i1];
        end
    end
%     StimulationFreqU=D.fsample/median(diff(stimulation));   %uncorrected stim frequency
    stimulation=stimulation(setdiff(1:length(stimulation),remove));
end

%remove isolated stimulations (short clusters)
dstimulation=diff(stimulation);
index = find(dstimulation>=1.4*Stim);
if ~isempty(index)
    index = unique([1 index+1 length(stimulation)+1]);
    dindex = diff(index);
    remove=[];
    for i1=2:length(dindex)
        if dindex(i1)<minStim
            remove=[remove sum(dindex(1:i1-1))+[1:dindex(i1)]];
        end
    end
%     StimulationFreqU=D.fsample/median(diff(stimulation));   %uncorrected stim frequency
    stimulation=stimulation(setdiff(1:length(stimulation),remove));
end



end
