function D=ImaGIN_InterpolationFilter(S)

% Correct artefacts in frequency space
% Ref: Ramirez & Baillet, 2011

% S.method = 'linear' or 'spline' or 'other'
if ~isfield(S, 'method') || isempty(S.method)
    S.method = 'pchip';
end

% try
    t=S.Fname;
% catch
%     t = spm_select(Inf, '\.mat$', 'Select data file');
% end
for i0=1:size(t,1)
    T=deblank(t(i0,:));
    D=spm_eeg_load(T);
    Time=time(D);
    
    
    F  = spm_figure('GetWin','Interactive');
    figure(F);clf
    
    
%     try
        EventType = S.EventType;
%     catch
%         EventType=spm_input('Events of interest', '+1', 'i',1);
%     end
%     try
        StartInterpolation = S.StartInterpolation;
%     catch
%         StartInterpolation=spm_input('Start of interpolation [sec]', '+1', 'r',-0.001);
%     end
%     try
        EndInterpolation = S.EndInterpolation;
%     catch
%         EndInterpolation=spm_input('End of interpolation [sec]', '+1', 'r',0.003);
%     end
    StartInterpolation=round(StartInterpolation*fsample(D));
    EndInterpolation=round(EndInterpolation*fsample(D));
    
    Bad=badchannels(D);
    Nchan=nchannels(D);
    Good=setdiff(1:Nchan,Bad);
    Nchan=length(Good);
    
    ev=events(D);
    
    Data=D(:,:);
    
    switch S.method
        
        case 'linear'
            
                %try to estimate artefact duration
            Ratio=1; %assumes artefact is Ratio times higher than the response
            startinterpolation=NaN*zeros(length(Good),length(ev));
            endinterpolation=NaN*zeros(length(Good),length(ev));
            indexint=[4*StartInterpolation:4*EndInterpolation];
            for i1=1:length(ev)
                if strcmp(EventType,strvcat(ev(i1).type))
                    ind=indsample(D,ev(i1).time);
                    ind1=[ind+4*StartInterpolation];
                    ind2=[ind+4*EndInterpolation];
                    Signal=Data(Good,ind1:ind2);
                    Signal=ImaGIN_normalisation(Signal,2,1:3*abs(StartInterpolation));
                    for i2=1:size(Signal,1)
                        Threshold=Ratio*mean(abs(Signal(i2,end-abs(EndInterpolation):end)));
                        tmp=find(abs(Signal(i2,:))>Threshold);
                        if ~isempty(tmp)
                            tmp2=min(find(tmp>=4*abs(StartInterpolation)&tmp<=4*abs(StartInterpolation)+2));
                            if ~isempty(tmp2)
                                tmpstart=tmp2;
                                tmpend=tmp2;
                                ok=1;
                                if tmp2==1
                                    ok=0;
                                end
                                while ok
                                    if isempty(intersect(tmp,tmp(tmpstart)-1))
                                        ok=0;
                                    else
                                        tmpstart=tmpstart-1;
                                    end
                                end
                                ok=1;
                                while ok
                                    if isempty(intersect(tmp,tmp(tmpend)+1))
                                        ok=0;
                                    else
                                        tmpend=tmpend+1;
                                    end
                                end
                                startinterpolation(i1,i2)=tmp(tmpstart);
                                endinterpolation(i1,i2)=tmp(tmpend);
                            end
                        end
                    end
                end
            end
            startinterpolation(startinterpolation==0)=NaN;
            endinterpolation(endinterpolation==0)=NaN;
            StartInterpolationEstimate=indexint(median(median(startinterpolation,2,'omitnan'),'omitnan'))-1;
            EndInterpolationEstimate=indexint(median(median(endinterpolation,2,'omitnan'),'omitnan'))+1;
            if StartInterpolationEstimate>StartInterpolation
                StartInterpolation=StartInterpolationEstimate;
            end
            if EndInterpolationEstimate<EndInterpolation
                EndInterpolation=EndInterpolationEstimate;
            end
            
            
            %Do the interpolation
        for i1=1:length(ev)
                if mod(i1,200)==0
                    disp([i1 length(ev)])
                end
                %         if ~isempty(intersect(EventType,ev(i1).value))
                if strcmp(EventType,strvcat(ev(i1).type))
                    ind=indsample(D,ev(i1).time);
                    tmp1=D(:,ind+[StartInterpolation EndInterpolation]);
                    tmp2=zeros(size(Data,1),length([StartInterpolation:EndInterpolation]));
                    for i2=1:size(Data,1)
                        tmp2(i2,:)=interp1([StartInterpolation EndInterpolation],tmp1(i2,:),[StartInterpolation:EndInterpolation]);
                    end
                    Data(:,ind+[StartInterpolation:EndInterpolation])=tmp2;
                end
            end
            
            
        case 'spline'
            
            %try to estimate artefact duration
            Ratio=5; %assumes artefact is Ratio times higher than the baseline (z-score)
            startinterpolation=NaN*zeros(length(Good),length(ev));
            endinterpolation=NaN*zeros(length(Good),length(ev));
            indexint=[4*StartInterpolation:4*EndInterpolation];
            GoodEvent=[];
            for i1=1:length(ev)
                if strcmp(EventType,strvcat(ev(i1).type))
                    GoodEvent=[GoodEvent i1];
                    ind=indsample(D,ev(i1).time);
                    ind1=[ind+4*StartInterpolation];
                    ind2=[ind+4*EndInterpolation];
                    Signal=Data(Good,ind1:ind2);
                    Signal=ImaGIN_normalisation(Signal,2,1:3*abs(StartInterpolation));
                    for i2=1:size(Signal,1)
%                         Threshold=Ratio*mean(abs(Signal(i2,end-abs(EndInterpolation):end)));
                        Threshold=Ratio;
                        tmp=find(abs(Signal(i2,:))>Threshold);
                        if ~isempty(tmp)
                            tmp2=min(find(tmp>=4*abs(StartInterpolation)-1&tmp<=4*abs(StartInterpolation)+3));
                            tmp3=max(find(tmp>=4*abs(StartInterpolation)-1&tmp<=4*abs(StartInterpolation)+3));
                            if ~isempty(tmp2)
                                tmpstart=tmp2;
                                tmpend=tmp3;
                                ok=1;
                                if tmp2==1
                                    ok=0;
                                end
                                while ok
                                    if isempty(intersect(tmp,[tmp(tmpstart)-1 tmp(tmpstart)-2]))
                                        ok=0;
                                    else
                                        tmpstart=tmpstart-1;
                                    end
                                end
                                ok=1;
                                while ok
                                    if isempty(intersect(tmp,[tmp(tmpend)+1 tmp(tmpend)+2]))
                                        ok=0;
                                    else
                                        tmpend=tmpend+1;
                                    end
                                end
                                startinterpolation(i2,i1)=tmp(tmpstart);
                                endinterpolation(i2,i1)=tmp(tmpend);
                            end
                        end
                    end
                end
            end
            startinterpolation=startinterpolation(:,GoodEvent);
            endinterpolation=endinterpolation(:,GoodEvent);

%             StartInterpolationEstimate=indexint(median(median(startinterpolation,2,'omitnan'),'omitnan'))-1;
%             EndInterpolationEstimate=indexint(median(median(endinterpolation,2,'omitnan'),'omitnan'))+1;
%             if StartInterpolationEstimate>StartInterpolation
%                 StartInterpolation=StartInterpolationEstimate;
%             end
%             if EndInterpolationEstimate<EndInterpolation
%                 EndInterpolation=EndInterpolationEstimate;
%             end
%    
%             %Do the interpolation
%             ind1=[];
%             ind2=[];
%             for i1=1:length(ev)
%                 %         if ~isempty(intersect(EventType,ev(i1).value))
%                 if ~isempty(intersect(EventType,ev(i1).type))
%                     ind=indsample(D,ev(i1).time);
%                     ind1=[ind1 ind+StartInterpolation];
%                     ind2=[ind2 ind+EndInterpolation];
%                 end
%             end
%             ind1=sort(ind1);
%             ind2=sort(ind2);
%             Time=time(D);
%             Index=1:length(Time);
%             for i1=1:length(ind1)
%                 Index=setdiff(Index,ind1(i1):ind2(i1));
%             end
%             Time2=Time(Index);
%             
%             for i1=1:size(Data,1)
%                 Data(i1,:) = spline(Time2,Data(i1,Index),Time);
%             end
            
                    %fill the NaNs and limit the range
            StartInterpolationEstimate=startinterpolation;
            EndInterpolationEstimate=endinterpolation;
            for i1=1:size(startinterpolation,1)
                if isnan(floor(median(startinterpolation(i1,:),2,'omitnan')))
                    startinterpolation(i1,1)=find(indexint==StartInterpolation);
                end
                if isnan(ceil(median(endinterpolation(i1,:),2,'omitnan')))
                    endinterpolation(i1,1)=find(indexint==EndInterpolation);
                end
                startinterpolation(i1,isnan(startinterpolation(i1,:)))=floor(median(startinterpolation(i1,:),2,'omitnan'));
                endinterpolation(i1,isnan(endinterpolation(i1,:)))=ceil(median(endinterpolation(i1,:),2,'omitnan'));
                StartInterpolationEstimate(i1,:)=indexint(startinterpolation(i1,:))-1;
                StartInterpolationEstimate(i1,find(StartInterpolationEstimate(i1,:)<StartInterpolation))=StartInterpolation;
                EndInterpolationEstimate(i1,:)=indexint(endinterpolation(i1,:))+1;
                EndInterpolationEstimate(i1,find(EndInterpolationEstimate(i1,:)>EndInterpolation))=EndInterpolation;
            end
            StartInterpolation=StartInterpolationEstimate;
            EndInterpolation=EndInterpolationEstimate;
            
            
            %Do the interpolation with event and channel specific values
            for i2=1:length(Good)
                
                ind1=[];
                ind2=[];
                n=0;
                for i1=1:length(ev)
                    %         if ~isempty(intersect(EventType,ev(i1).value))
                    if ~isempty(intersect(EventType,ev(i1).type))
                        n=n+1;
                        ind=indsample(D,ev(i1).time);
                        ind1=[ind1 ind+StartInterpolation(i2,n)];
                        ind2=[ind2 ind+EndInterpolation(i2,n)];
                    end
                end
                ind1=sort(ind1);
                ind2=sort(ind2);
                Time=time(D);
                Index=1:length(Time);
                for i1=1:length(ind1)
                    Index=setdiff(Index,ind1(i1):ind2(i1));
                end
                Time2=Time(Index);
                
                Data(Good(i2),:) = spline(Time2,Data(Good(i2),Index),Time);
            end
            
    
            
            
        case 'pchip'
            
            %try to estimate artefact duration
            Ratio=5; %assumes artefact is Ratio times higher than the baseline (z-score)
            startinterpolation=NaN*zeros(length(Good),length(ev));
            endinterpolation=NaN*zeros(length(Good),length(ev));
            indexint=[4*StartInterpolation:4*EndInterpolation];
            GoodEvent=[];
            for i1=1:length(ev)
                disp(ev(i1))                
                if strcmp(EventType,strvcat(ev(i1).type))
                    GoodEvent=[GoodEvent i1];
                    ind=indsample(D,ev(i1).time);
                    ind1=[ind+4*StartInterpolation];
                    ind2=[ind+4*EndInterpolation];
                    Signal=Data(Good,ind1:ind2);
                    Signal=ImaGIN_normalisation(Signal,2,1:3*abs(StartInterpolation));
                    for i2=1:size(Signal,1)
%                         Threshold=Ratio*mean(abs(Signal(i2,end-abs(EndInterpolation):end)));
                        Threshold=Ratio;
                        tmp=find(abs(Signal(i2,:))>Threshold);
                        if ~isempty(tmp)
                            tmp2=min(find(tmp>=4*abs(StartInterpolation)-1&tmp<=4*abs(StartInterpolation)+3));
                            tmp3=max(find(tmp>=4*abs(StartInterpolation)-1&tmp<=4*abs(StartInterpolation)+3));
                            if ~isempty(tmp2)
                                tmpstart=tmp2;
                                tmpend=tmp3;
                                ok=1;
                                if tmp2==1
                                    ok=0;
                                end
                                while ok
                                    if isempty(intersect(tmp,[tmp(tmpstart)-1 tmp(tmpstart)-2]))
                                        ok=0;
                                    else
                                        tmpstart=tmpstart-1;
                                    end
                                end
                                ok=1;
                                while ok
                                    if isempty(intersect(tmp,[tmp(tmpend)+1 tmp(tmpend)+2]))
                                        ok=0;
                                    else
                                        tmpend=tmpend+1;
                                    end
                                end
                                startinterpolation(i2,i1)=tmp(tmpstart);
                                endinterpolation(i2,i1)=tmp(tmpend);
                            end
                        end
                    end
                end
            end
            startinterpolation=startinterpolation(:,GoodEvent);
            endinterpolation=endinterpolation(:,GoodEvent);
            
            
%             StartInterpolationEstimate=indexint(median(median(startinterpolation,2,'omitnan'),'omitnan'))-1;
%             EndInterpolationEstimate=indexint(median(median(endinterpolation,2,'omitnan'),'omitnan'))+1;
%             if StartInterpolationEstimate>StartInterpolation
%                 StartInterpolation=StartInterpolationEstimate;
%             end
%             if EndInterpolationEstimate<EndInterpolation
%                 EndInterpolation=EndInterpolationEstimate;
%             end
%    
%             %Do the interpolation
%             ind1=[];
%             ind2=[];
%             for i1=1:length(ev)
%                 %         if ~isempty(intersect(EventType,ev(i1).value))
%                 if ~isempty(intersect(EventType,ev(i1).type))
%                     ind=indsample(D,ev(i1).time);
%                     ind1=[ind1 ind+StartInterpolation];
%                     ind2=[ind2 ind+EndInterpolation];
%                 end
%             end
%             ind1=sort(ind1);
%             ind2=sort(ind2);
%             Time=time(D);
%             Index=1:length(Time);
%             for i1=1:length(ind1)
%                 Index=setdiff(Index,ind1(i1):ind2(i1));
%             end
%             Time2=Time(Index);
%             
%             for i1=1:size(Data,1)
%                 Data(i1,:) = pchip(Time2,Data(i1,Index),Time);
%             end


            %fill the NaNs and limit the range
            StartInterpolationEstimate=startinterpolation;
            EndInterpolationEstimate=endinterpolation;
            for i1=1:size(startinterpolation,1)
                if isnan(floor(median(startinterpolation(i1,:),2,'omitnan')))
                    startinterpolation(i1,1)=find(indexint==StartInterpolation);
                end
                if isnan(ceil(median(endinterpolation(i1,:),2,'omitnan')))
                    endinterpolation(i1,1)=find(indexint==EndInterpolation);
                end
                startinterpolation(i1,isnan(startinterpolation(i1,:)))=floor(median(startinterpolation(i1,:),2,'omitnan'));
                endinterpolation(i1,isnan(endinterpolation(i1,:)))=ceil(median(endinterpolation(i1,:),2,'omitnan'));
                StartInterpolationEstimate(i1,:)=indexint(startinterpolation(i1,:))-1;
                StartInterpolationEstimate(i1,find(StartInterpolationEstimate(i1,:)<StartInterpolation))=StartInterpolation;
                EndInterpolationEstimate(i1,:)=indexint(endinterpolation(i1,:))+1;
                EndInterpolationEstimate(i1,find(EndInterpolationEstimate(i1,:)>EndInterpolation))=EndInterpolation;
            end
            StartInterpolation=StartInterpolationEstimate;
            EndInterpolation=EndInterpolationEstimate;
            
            
            %Do the interpolation with event and channel specific values
            for i2=1:length(Good)
                
                ind1=[];
                ind2=[];
                n=0;
                for i1=1:length(ev)
                    %         if ~isempty(intersect(EventType,ev(i1).value))
                    if strcmp(EventType,ev(i1).type)
                        n=n+1;
                        ind=indsample(D,ev(i1).time);
                        ind1=[ind1 ind+StartInterpolation(i2,n)];
                        ind2=[ind2 ind+EndInterpolation(i2,n)];
                    end
                end
                ind1=sort(ind1);
                ind2=sort(ind2);
                Time=time(D);
                Index=1:length(Time);
                for i1=1:length(ind1)
                    Index=setdiff(Index,ind1(i1):ind2(i1));
                end
                Time2=Time(Index);
                
                Data(Good(i2),:) = pchip(Time2,Data(Good(i2),Index),Time);
            end
            
        case 'svd'
            
            %Define basis function
            PeakMin=-0.001;
            PeakMax=0.003;      %Could be defined as argument
            Artefact=[];
            for i1=1:length(ev)
                if strcmp(EventType,strvcat(ev(i1).type))
                    T2=indsample(D,S.StartInterpolation+ev(i1).time);
                    T3=indsample(D,S.EndInterpolation+ev(i1).time);
                    T4=indsample(D,PeakMin+ev(i1).time);
                    T5=indsample(D,PeakMax+ev(i1).time);
                    T6=indsample(D,ev(i1).time);
                    Signe=sign(mean(D(:,(T6-1):(T6+1)),2)-mean(D(:,T2:T3),2));
                    if ~isnan(T3)
                        Artefact=cat(1,Artefact,(Signe(Good)*ones(1,T3-T2+1)).*D(Good,T2:T3));
                    end
                end
            end
            [u,s,v]=svd(Artefact,'econ');
            V=abs(v-ones(size(v,2),1)*v(1,:));
            index=[T4:T5]-T2+1;
            
%             %Selection based on vaiance
%             tmp=cumsum(diag(s))/sum(diag(s));
%             Select=1:max(find(tmp<0.99));

            %Selection based on the max
            [Y,I] = max(V,[],1);
            Select=find(I>=min(index)&I<=max(index));
            
%             %Selection based on the ratio
%             Ratio=sum(V(index,:))./sum(V);
%             Select=find(Ratio>length(index)/(size(V,1)-1));

            Basis=[v(:,Select)];
            
            
%             %Do local correction to components
%             for i1=1:size(v,2)
%                 k = sort([local_max(v(:,i1))' local_max(-v(:,i1))']);
%                 indexk=intersect(index,k);
%                 if ~isempty(indexk)
% %                     v([(index(1)-1):(index(end)+1)],i1)=interp1([index(1)-1 index(end)+1],v([index(1)-1 index(end)+1],i1),[(index(1)-1):(index(end)+1)]);
%                     v(:,i1)=spline([1:index(1)-1 index(end)+1:size(v,1)],v([1:index(1)-1 index(end)+1:size(v,1)],i1),[1:size(v,1)]);
%                 end
%             end
%             Basis=v;
           
            
            %Apply correction
            n=0;
            for i1=1:length(ev)
                if strcmp(EventType,strvcat(ev(i1).type))
                    T2=indsample(D,S.StartInterpolation+ev(i1).time);
                    T3=indsample(D,S.EndInterpolation+ev(i1).time);
                    if ~isnan(T3)
                        Correction=D(:,T2:T3)*Basis*Basis';
                        Data(:,T2:T3)=D(:,T2:T3)-Correction;
                        Bias=(Data(:,T2)+Data(:,T3)-D(:,T2-1)-D(:,T3+1))./2;
                        Data(:,T2:T3)=Data(:,T2:T3)-Bias;
%                     T2=indsample(D,S.StartInterpolation+ev(i1).time);
%                     T3=indsample(D,S.EndInterpolation+ev(i1).time);
%                     T4=indsample(D,PeakMin+ev(i1).time);
%                     T5=indsample(D,PeakMax+ev(i1).time);
%                     T6=indsample(D,ev(i1).time);
%                     Signe=sign(mean(D(:,(T6-1):(T6+1)),2)-mean(D(:,T2:T3),2));
%                         n=n+1;
%                         indexchan=(n-1)*length(Good)+[1:length(Good)];
%                         Data(Good,T2:T3)=(Signe(Good)*ones(1,T3-T2+1)).*u(indexchan,:)*s*Basis';
                    end
                end
            end

%             %Define basis function
%             PeakMin=-0.001;
%             PeakMax=0.004;      %Could be defined as argument
%             Artefact=[];
%             %resample data at 2 kHz
% %             clear SS
% %             SS.D=S.Fname;
% %             SS.fsample_new=4e3;
% %             Ds = spm_eeg_downsample(SS);
% %             Data=D(:,:);
%             timeS=min(time(D)):(1/2000):max(time(D));
%             Ds=zeros(nchannels(D),length(timeS));
%             for i1=1:nchannels(D)
%                 Ds(i1,:)=interp1(time(D),D(i1,:),timeS);
%             end
% 
%             for i1=1:length(ev)
%                 if strcmp(EventType,strvcat(ev(i1).type))
%                     T2=min(find(abs(timeS-(S.StartInterpolation+ev(i1).time))==min(abs(timeS-(S.StartInterpolation+ev(i1).time)))));
%                     T3=min(find(abs(timeS-(S.EndInterpolation+ev(i1).time))==min(abs(timeS-(S.EndInterpolation+ev(i1).time)))));
%                     T4=min(find(abs(timeS-(PeakMin+ev(i1).time))==min(abs(timeS-(PeakMin+ev(i1).time)))));
%                     T5=min(find(abs(timeS-(PeakMax+ev(i1).time))==min(abs(timeS-(PeakMax+ev(i1).time)))));
%                     T6=min(find(abs(timeS)==min(abs(timeS))));
%                     Signe=sign(mean(Ds(:,(T6-1):(T6+1)),2)-mean(Ds(:,T2:T3),2));
%                     if ~isnan(T3)
%                         Artefact=cat(1,Artefact,(Signe(Good)*ones(1,T3-T2+1)).*Ds(Good,T2:T3));
%                     end
%                 end
%             end
%             [u,s,v]=svd(Artefact,'econ');
%             V=abs(v-ones(size(v,2),1)*v(1,:));
%             index=[T4:T5]-T2+1;
%             
% %             %Selection based on vaiance
% %             tmp=cumsum(diag(s))/sum(diag(s));
% %             Select=1:max(find(tmp<0.99));
%             
%             %Selection based on the max
%             [Y,I] = max(V,[],1);
%             Select=find(I>=min(index)&I<=max(index));
%             
% %             %Selection based on the ratio
% %             Ratio=sum(V(index,:))./sum(V);
% %             Select=find(Ratio>length(index)/(size(V,1)-1));
% 
%             Basis=[v(:,Select)];
%             
%             %Apply correction
%             Datas=Ds(:,:);
%             for i1=1:length(ev)
%                 if strcmp(EventType,strvcat(ev(i1).type))
%                     T2=min(find(abs(timeS-(S.StartInterpolation+ev(i1).time))==min(abs(timeS-(S.StartInterpolation+ev(i1).time)))));
%                     T3=min(find(abs(timeS-(S.EndInterpolation+ev(i1).time))==min(abs(timeS-(S.EndInterpolation+ev(i1).time)))));
%                     if ~isnan(T3)
%                         Correction=Ds(:,T2:T3)*Basis*Basis';
%                         Datas(:,T2:T3)=Ds(:,T2:T3)-Correction;
%                         Bias=(Datas(:,T2)+Datas(:,T3)-Ds(:,T2-1)-Ds(:,T3+1))./2;
%                         Datas(:,T2:T3)=Datas(:,T2:T3)-Bias;
%                     end
%                 end
%             end
%             
%             %Reinterpolate
%             Data=D(:,:);
%             for i1=1:nchannels(D)
%                 Data(i1,:)=interp1(timeS,Datas(i1,:),time(D));
%             end
%             
%             
% %     figure
% %     for i2=1:nchannels(Ds)
% %         plot([Ds(i2,:);Datas(i2,:)]')
% %         title(num2str(i2))
% %         pause
% %     end

            
        case 'other' 
            
            NCompo=2;
            CCthresh=0.5;
            ind12=[];
            ind22=[];
            for i1=1:length(ev)
                %         if ~isempty(intersect(EventType,ev(i1).value))
                if strcmp(EventType,strvcat(ev(i1).type))
                    ind=indsample(D,ev(i1).time);
                    ind12=[ind12 ind+StartInterpolation*3];
                    ind22=[ind22 ind+EndInterpolation*3];
                end
            end
            ind12=sort(ind12);
            ind22=sort(ind22);
            Index2=[];
            for i1=1:length(ind12)
                Index2=[Index2 ind12(i1):ind22(i1)];
            end
            [u s v]=svd(Data(:,Index2)',0);
            
            f=abs(fft(u));
            CC=zeros(1,size(f,2));
            for i1=1:size(f,2)
                cc=corrcoef(f(:,1),f(:,i1));
                CC(i1)=cc(2);
            end
            CompoArtefact=find(CC>CCthresh);
            %     s=diag(s);
            %     s(CompoArtefact)=0;
            %     s=diag(s);
            
            
            U=u(:,CompoArtefact);
            windowSize = round(1*size(U,1)/length(ev));
            for i1=1:length(CompoArtefact)
                U2(:,i1)=2*U(:,i1)-filter(ones(1,windowSize)/windowSize,1,U(:,i1))-flipud(filter(ones(1,windowSize)/windowSize,1,flipud(U(:,i1))));
            end
            Artefact=ImaGIN_normalisation(sum(abs(U2),2),1);
            windowSize = round(1*length(Artefact)/length(ev));
            Artefact=2*Artefact-filter(ones(1,windowSize)/windowSize,1,Artefact)-flipud(filter(ones(1,windowSize)/windowSize,1,flipud(Artefact)));
            
            Index3=setdiff(1:length(Artefact),setdiff(find(Artefact>0.5),[1 length(Artefact)]));
            Time3=Time(Index2);
            Time4=Time3(Index3);
            %     for i1=1:NCompo
            for i1=CompoArtefact
                %         u(:,i1) = spline(Time4,u(Index3,i1),Time3);
                u(:,i1) = interp1(Time4,u(Index3,i1),Time3,'linear');
            end
            Data(:,Index2)=transpose(u*s*v');
            
            
    end
    
    
    %Save as a newfile
    D=clone(D,['i' D.fname], [D.nchannels D.nsamples D.ntrials]);
    D(:,:,:)=Data;
    save(D);
end
