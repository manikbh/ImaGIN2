dirFTRACT = '/gin/data/database/03-preprocessed/FTRACT/';

pCODES = dir(fullfile(dirFTRACT,'0*'));
pCODES = {pCODES.name};

for i = 1 : length(pCODES)
    dirSEEG = fullfile(dirFTRACT,pCODES{i}, 'SEEG');
    if exist(dirSEEG,'dir')
        interICTAL = dir(fullfile(dirSEEG,'interICTAL*'));
        interICTAL  = {interICTAL.name};
        if ~isempty(interICTAL)
            dirSPK = fullfile(dirSEEG,interICTAL{1}); 
            idx_sc = strfind(dirSPK,'_');
            if ~isempty(idx_sc)
                seeg_date = dirSPK(idx_sc(end)+1:end);
            end
            SPK_files = dir(fullfile(dirSPK,'SPK_*.mat'));
            SPK_files = {SPK_files.name};            
            if ~isempty(SPK_files)  
                fprintf(fid,'%s\n', strcat(pCODES{i},'_',seeg_date));
                cumulRate = 0;
                nf = 0;
                for j = 1:length(SPK_files)
                    SPK_file = fullfile(dirSPK,SPK_files{j}); 
                    load(SPK_file);
                    spikeRateFile = spkFile.spikeRate;
                    newLabels = ImaGIN_channels_standard(spkFile.spikeRate.labels,'all_upper', 'bipolar');
                    spikeRateFile.labels = newLabels;
                    spikeRate     = table2cell(spikeRateFile);
                    if exist(fullfile(dirSPK,'spikeRate.mat'),'file')
                        delete(fullfile(dirSPK,'spikeRate.mat'));
                    end
                    savedFile = fullfile(dirSPK,['spikeRate_', num2str(j),'.mat']);
                    save(savedFile,'spikeRate');
                    try
                        cumulRate     = cumulRate + spikeRateFile.spk_Rate;
                        nf            = nf + 1;
                    catch
                        disp(['WARNING: ' pCODES{i} ' recorded several dates !'])
                    end
                end
                if nf == numel(SPK_files)
                    meanSPK = cumulRate / nf;
                    spikeRateFile.spk_Rate = meanSPK;
                    savedFileSPK = fullfile(dirSPK,['spikeRate_','global','.mat']);
                    globalSPK    = table2cell(spikeRateFile);
                    save(savedFileSPK,'globalSPK');
                    fprintf([savedFileSPK, '  is created \n']);
                    
                end
            end
        end       
    end  
end

