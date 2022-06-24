function prepare_ImaGIN_ArtefactCorrectionModel(tms, time_break, rc1, rc2, nb, mode, channels, adjustement, screenshot, estimate, FileIn, FileOut, DirScreenshots)

% FileIn: path linking to the MEEG file to correct
% DirOut: output path
% tms: total stimulation time (s)
% rc1: minimum value of RC (s) (default:1e-5)
% rc2: maximum value of RC (s) (default:1e-3)
% nb: number of artifact shape generated (default:500)
% mode: stimulation mode (monophasic or biphasic)

S.Fname=FileIn;

if ischar(tms)
    tms=str2num(tms);
end
S.tms=tms;

if ischar(time_break)
    time_break=str2num(time_break);
end
S.break=time_break;

if ischar(rc1)
    rc1=str2num(rc1);
end
S.rc1=rc1;

if ischar(rc2)
    rc2=str2num(rc2);
end
S.rc2=rc2;

if ischar(nb)
    nb=str2num(nb);
end
S.nb=nb;

if ischar(channels)
    channels = str2num(channels);
    if channels == 0
        channels = [];
    end
end
S.channels = channels;

if strcmp(screenshot,'False')
    screenshot=false;
else if strcmp(screenshot,'True')
        screenshot=true;
    end
end
S.screenshot = screenshot;

if strcmp(estimate,'False')
    estimate=false;
else if strcmp(estimate,'True')
        estimate=true;
    end
end

if ischar(adjustement)
    adjustement=str2num(adjustement);
end
S.adjustement=adjustement;

S.artefactEstimation = estimate;

S.DirScreenshots = DirScreenshots;

S.mode=mode;
S.FileOut=FileOut;
D=ImaGIN_ArtefactCorrectionModel(S);

set_final_status('OK')

close all

end
