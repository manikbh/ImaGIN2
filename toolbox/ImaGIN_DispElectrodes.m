function ImaGIN_DispElectrodes(S)
% Display electrodes positions.
%
% INPUTS:
%    - S: Structure with optional fields
%      |- Fname : Full path to the .mat/.dat file to display
%      |- P     : Full path to the T1 overlay to use in the display

% -=============================================================================
% This function is part of the ImaGIN software: 
% https://f-tract.eu/
%
% This software is distributed under the terms of the GNU General Public License
% as published by the Free Software Foundation. Further details on the GPLv3
% license can be found at http://www.gnu.org/copyleft/gpl.html.
%
% FOR RESEARCH PURPOSES ONLY. THE SOFTWARE IS PROVIDED "AS IS," AND THE AUTHORS
% DO NOT ASSUME ANY LIABILITY OR RESPONSIBILITY FOR ITS USE IN ANY CONTEXT.
%
% Copyright (c) 2000-2018 Inserm U1216
% =============================================================================-
%
% Authors: Olivier David, Francois Tadel


% ===== PARSE INPUTS =====
% Get file to display
if (nargin >= 1) && isfield(S, 'Fname') && ~isempty(S.Fname)
    t = S.Fname;
else
    t = spm_select(1, '\.mat$', 'Select SEEG data file');
    if isempty(t)
        return;
    end
end
% T1 overlay
if (nargin >= 1) && isfield(S, 'Fname') && ~isempty(S.Fname)
    P = S.P;
else
    P = [];
end


% ===== LOAD SEEG POSITIONS =====
% Load file
D = spm_eeg_load(t);
% Get electrodes names and labels
Sensors = sensors(D,'EEG');
% Check for missing information
if isempty(Sensors) || isempty(Sensors.elecpos)
    error('No electrode positions in this file.');
end
% Exclude NaN and (0,0,0) positions 
iGood = find(~any(isnan(Sensors.elecpos),2) & ~all(Sensors.elecpos == 0,2));
if isempty(iGood)
    error('The electrodes positions have not been defined in this file.');
end
% Get positions and labels of the good channels
chLoc   = Sensors.elecpos(iGood,:);
chNames = Sensors.label(iGood);


% ===== GET VOLUME OVERLAY =====
if isempty(P)
    % Template or individual subject?
    tmp = spm_input('T1 overlay ', '+1', 'Template|Subject');
    % Template
    if strcmp(tmp,'Template')
        % Ask for type of atlas to use
        Species = spm_input('Species ','+1','Human|Rat|Mouse');
        % Select corresponding atlas file
        switch Species
            case 'Human'
                P = fullfile(spm('dir'),'canonical','avg152T1.nii');
            case 'Rat'
                P = fullfile(spm_str_manip(spm('dir'),'h'),'ratlas5','template','template_T1.img');
            case 'Mouse' 
                P = fullfile(spm_str_manip(spm('dir'),'h'),'mousatlas5','template','template_T1.img');
        end
    % Individual brain
    else
        P = spm_select(1, 'image', 'Select image for rendering on');
    end
end

% ===== DISPLAY ELECTRODES =====
Nchannels = size(chLoc,1);
% Prepare input structures
sdip.n_seeds = 1;
sdip.n_dip   = Nchannels;
sdip.Mtb     = 1;
sdip.j{1}    = zeros(3 * Nchannels, 1);
sdip.loc{1}  = chLoc';
sdip.Names   = chNames;
sdip.fig     = spm_figure('GetWin','Graphics');
% Render electrodes
ImaGIN_DrawElectrodes('Init', sdip, P);

end



%% ================================================================================================
function varargout = ImaGIN_DrawElectrodes(action,varargin)
% ImaGIN_DrawElectrodes
%
% Function to display the dipoles as obtained from the optim routine.
%
% Use it with arguments or not:
% - ImaGIN_DrawElectrodes('Init')
%        The routine asks for the dipoles file and image to display
% - ImaGIN_DrawElectrodes('Init',sdip)
%        The routine will use the avg152T1 canonical image 
% - ImaGIN_DrawElectrodes('Init',sdip,P)
%        The routines dispays the dipoles on image P.
%
% If multiple seeds have been used, you can select the seeds to display 
% by pressing their index. 
% Given that the sources could have different locations, the slices
% displayed will be the 3D view at the *average* or *mean* locations of
% selected sources.
% If more than 1 dipole was fitted at a time, then selection of source 1 
% to N is possible through the pull-down selector.
%
% The location of the source/cut is displayed in mm and voxel, as well as
% the underlying image intensity at that location.
% The cross hair position can be hidden by clicking on its button.
%
% Nota_1: If the cross hair is manually moved by clicking in the image or
%       changing its coordinates, the dipole displayed will NOT be at
%       the right displayed location. That's something that needs to be improved...
%
% Nota_2: Some seeds may have not converged within the limits fixed,
%       these dipoles are not displayed...
%
% Fields needed in sdip structure to plot on an image:
%       + n_seeds: nr of seeds set used, i.e. nr of solutions calculated
%       + n_dip: nr of fitted dipoles on the EEG time series
%       + loc: location of fitted dipoles, cell{1,n_seeds}(3 x n_dip)
%               remember that loc is fixed over the time window.
%       + j: sources amplitude over the time window, 
%            cell{1,n_seeds}(3*n_dip x Ntimebins)
%       + Mtb: index of maximum power in EEG time series used
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Christophe Phillips,
% $Id$

global st

Fig = spm_figure('FindWin');
colors = strvcat('y','b','g','r','c','m');              % 6 possible colors
marker = strvcat('o','x','+','*','s','d','v','p','h');  % 9 possible markers
Ncolors = length(colors);
Nmarker = length(marker);

if nargin == 0, action = 'Init'; end;

switch lower(action),
%________________________________________________________________________
case 'init'
%------------------------------------------------------------------------
% FORMAT ImaGIN_DrawElectrodes('Init',sdip,P)
% Initialise the variables with GUI
% e.g. ImaGIN_DrawElectrodes('Init')
%------------------------------------------------------------------------
spm('Clear')

if nargin<2
    load(spm_select(1,'^.*S.*dip.*\.mat$','Select dipole file'));
    if ~exist('sdip') && exist('result')
        sdip = result;
    end
else
    sdip = varargin{1};
end

% if the exit flag is not in the structure, assume everything went ok.
if ~isfield(sdip,'exitflag')
    sdip.exitflag = ones(1,sdip.n_seeds);
end

Pcanonical = fullfile(spm('dir'),'canonical','avg152T1.nii');

if nargin<3
	if ~isstruct(st)
        P = Pcanonical;
	elseif isempty(st.vols{1})
        P = Pcanonical;
	end
else
    P = varargin{2};
end
if ischar(P), P = spm_vol(P); end;
spm_orthviews('Reset');
spm_orthviews('Image', P, [0.0 0.45 1 0.55]);
spm_orthviews('Interp', 0); %OD
spm_orthviews('MaxBB');
st.callback = 'spm_image(''shopos'');';


WS = spm('WinScale');

% Build GUI
%=============================
% Location:
%-----------
uicontrol(Fig,'Style','Frame','Position',[60 25 200 325].*WS,'DeleteFcn','spm_image(''reset'');');
uicontrol(Fig,'Style','Frame','Position',[70 250 180 90].*WS);
uicontrol(Fig,'Style','Text', 'Position',[75 320 170 016].*WS,'String','Crosshair Position');
uicontrol(Fig,'Style','PushButton', 'Position',[75 316 170 006].*WS,...
	'Callback','spm_orthviews(''Reposition'',[0 0 0]);','ToolTipString','move crosshairs to origin');
uicontrol(Fig,'Style','Text', 'Position',[75 295 35 020].*WS,'String','mm:');
uicontrol(Fig,'Style','Text', 'Position',[75 275 35 020].*WS,'String','vx:');
uicontrol(Fig,'Style','Text', 'Position',[75 255 75 020].*WS,'String','Img Intens.:');

st.mp = uicontrol(Fig,'Style','edit', 'Position',[110 295 135 020].*WS,'String','','Callback','spm_image(''setposmm'')','ToolTipString','move crosshairs to mm coordinates','Tag','spm_image:mm');
st.vp = uicontrol(Fig,'Style','edit', 'Position',[110 275 135 020].*WS,'String','','Callback','spm_image(''setposvx'')','ToolTipString','move crosshairs to voxel coordinates','Tag','spm_image:vx');
st.in = uicontrol(Fig,'Style','Text', 'Position',[150 255  85 020].*WS,'String','','Tag','spm_image:intensity');

c = 'if get(gco,''Value'')==1, spm_orthviews(''Xhairs'',''off''), else, spm_orthviews(''Xhairs'',''on''); end;';
uicontrol(Fig,'Style','togglebutton','Position',[95 220 125 20].*WS,...
	'String','Hide Crosshairs','Callback',c,'ToolTipString','show/hide crosshairs');

% Dipoles/seeds selection:
%--------------------------
uicontrol(Fig,'Style','Frame','Position',[300 25 180 325].*WS);
uicontrol(Fig,'Style','text','String','Select channel:','HorizontalAlignment','left', ...
        'Position',[310 320 90 20].*WS);
txt_box = cell(sdip.n_dip,1);
for ii=1:sdip.n_dip, txt_box{ii} = sdip.Names{ii}; end
txt_box{sdip.n_dip+1} = 'all';
sdip.hdl.hdip = uicontrol(Fig,'Style','popup','String',txt_box, ...
        'Position',[400 323 75 20].*WS, ...
        'Callback','ImaGIN_DrawElectrodes(''ChgDip'')');
sdip.hdl.hdipplus = uicontrol(Fig,'Style','pushbutton','String','+', ...
        'Position',[450 293 20 20].*WS, ...
        'Callback','ImaGIN_DrawElectrodes(''ChgDip'')');
sdip.hdl.hdipminus = uicontrol(Fig,'Style','pushbutton','String','-', ...
        'Position',[400 293 20 20].*WS, ...
        'Callback','ImaGIN_DrawElectrodes(''ChgDip'')','Visible','off');
        
% Dipoles orientation and strength:
%-----------------------------------
uicontrol(Fig,'Style','Frame','Position',[70 120 180 90].*WS);
uicontrol(Fig,'Style','Text', 'Position',[75 190 170 016].*WS, ...
        'String','Dipole orientation & strength');
uicontrol(Fig,'Style','Text', 'Position',[75 165 65 020].*WS,'String','vn_x_y_z:');
uicontrol(Fig,'Style','Text', 'Position',[75 145 75 020].*WS,'String','theta, phi:');
uicontrol(Fig,'Style','Text', 'Position',[75 125 75 020].*WS,'String','Dip intens.:');
   
sdip.hdl.hor1 = uicontrol(Fig,'Style','Text', 'Position',[140 165  105 020].*WS,'String','a');
sdip.hdl.hor2 = uicontrol(Fig,'Style','Text', 'Position',[150 145  85 020].*WS,'String','b');
sdip.hdl.int  = uicontrol(Fig,'Style','Text', 'Position',[150 125  85 020].*WS,'String','c');
    
st.vols{1}.sdip = sdip;

% First plot = first channel
ImaGIN_DrawElectrodes('DrawDip',1,1)
    
%________________________________________________________________________
case 'drawdip'
%------------------------------------------------------------------------
% FORMAT ImaGIN_DrawElectrodes('DrawDip',i_seed,i_dip,sdip)
% e.g. ImaGIN_DrawElectrodes('DrawDip',1,1,sdip)
% e.g. ImaGIN_DrawElectrodes('DrawDip',[1:5],1,sdip)

% when displaying a set of 'close' dipoles
% this defines the limit between 'close' and 'far'
% lim_cl = 10; %mm
% No more limit. All source displayed as projected on mean 3D cut.

if nargin<2
    error('At least i_seed')
end
i_seed = varargin{1};
if nargin<3
    i_dip = 1;
else
    i_dip = varargin{2}; 
end

if nargin<4
    if isfield(st.vols{1},'sdip')
        sdip = st.vols{1}.sdip;
    else
        error('I can''t find sdip structure');
    end
else
    sdip = varargin{3};
    st.vols{1}.sdip = sdip;
end

if any(i_seed>sdip.n_seeds) | i_dip>(sdip.n_dip+1)
    error('Wrong i_seed or i_dip index in ImaGIN_DrawElectrodes');
end

% Note if i_dip==(sdip.n_dip+1) all dipoles are displayed simultaneously
if i_dip==(sdip.n_dip+1)
    i_dip=1:sdip.n_dip;
end

% if seed indexes passed is wrong (no convergence) remove the wrong ones
i_seed(find(sdip.exitflag(i_seed)~=1)) = [];
if isempty(i_seed)
    error('You passed the wrong seed indexes...')
end
if size(i_seed,2)==1, i_seed=i_seed'; end

% Display business
%-----------------
loc_mm = sdip.loc{i_seed(1)}(:,i_dip);
if length(i_seed)>1
    for ii = i_seed(2:end)
        loc_mm = loc_mm + sdip.loc{ii}(:,i_dip);
    end
    loc_mm = loc_mm/length(i_seed);
end
if length(i_dip)>1
    loc_mm = mean(loc_mm,2);
end

% Place the underlying image at right cuts
spm_orthviews('Reposition',loc_mm);

if length(i_dip)>1
    tabl_seed_dip = [kron(ones(length(i_dip),1),i_seed') ...
                        kron(i_dip',ones(length(i_seed),1))];
else
    tabl_seed_dip = [i_seed' ones(length(i_seed),1)*i_dip];
end

% Scaling, according to all dipoles in the selected seed sets.
% The time displayed is the one corresponding to the maximum EEG power !
Mn_j = -1;
l3 = -2:0;
for ii = 1:length(i_seed)
    for jj = 1:sdip.n_dip
        Mn_j = max([Mn_j sqrt(sum(sdip.j{ii}(jj*3+l3,sdip.Mtb).^2))]);
    end
end
if Mn_j==0
    Mn_j=1;
end
    
st.vols{1}.sdip.tabl_seed_dip = tabl_seed_dip;

% Display all dipoles, the 1st one + the ones close enough.
% Run through the 6 colors and 9 markers to differentiate the dipoles.
% NOTA: 2 dipoles coming from the same set will have same colour/marker
ind = 1 ;
dip_h = zeros(6,size(tabl_seed_dip,1),1);
    % each dipole displayed has 6 handles:
    %   2 per view (2*3): one for the line, the other for the circle
js_m = zeros(3,1);

% Deal with case of multiple i_seed and i_dip displayed.
% make sure dipole from same i_seed have same colour but different marker.

pi_dip = find(diff(tabl_seed_dip(:,2)));
if isempty(pi_dip)
    % i.e. only one dip displayed per seed, use old fashion
    for ii=1:size(tabl_seed_dip,1)
        if ii>1
            if tabl_seed_dip(ii,1)~=tabl_seed_dip(ii-1,1)
                ind = ind+1;
            end
        end
        ic = mod(ind-1,Ncolors)+1;
        im = fix(ind/Ncolors)+1;

        loc_pl = sdip.loc{tabl_seed_dip(ii,1)}(:,tabl_seed_dip(ii,2));

        js = sdip.j{tabl_seed_dip(ii,1)}(tabl_seed_dip(ii,2)*3+l3,sdip.Mtb);
        dip_h(:,ii) = add1dip(loc_pl,js/Mn_j*20,marker(im),colors(ic),st.vols{1}.ax,Fig,st.bb);
        js_m = js_m+js;
    end
else
    for ii=1:pi_dip(1)
        if ii>1
            if tabl_seed_dip(ii,1)~=tabl_seed_dip(ii-1,1)
                ind = ind+1;
            end
        end
        for jj=1:sdip.n_dip
            im = mod(jj-1,Nmarker)+1;
            ic = mod(jj-1,Ncolors)+1;
            
            loc_pl = sdip.loc{tabl_seed_dip(ii,1)}(:,jj);
            js = sdip.j{tabl_seed_dip(ii,1)}(jj*3+l3,sdip.Mtb);
            js_m = js_m+js;
            dip_h(:,ii+(jj-1)*pi_dip(1)) = ...
                add1dip(loc_pl,js/Mn_j*20,marker(im),colors(ic), ...
                        st.vols{1}.ax,Fig,st.bb);
        end
    end
end
st.vols{1}.sdip.ax = dip_h;

% Display dipoles orientation and strength
js_m = js_m/size(tabl_seed_dip,1);
[th,phi,Ijs_m] = cart2sph(js_m(1),js_m(2),js_m(3));
if Ijs_m==0
    Ijs_m=1;
end
Njs_m = round(js_m'/Ijs_m*100)/100;
Angle = round([th phi]*1800/pi)/10;

set(sdip.hdl.hor1,'String',[num2str(Njs_m(1)),' ',num2str(Njs_m(2)), ...
        ' ',num2str(Njs_m(3))]);
set(sdip.hdl.hor2,'String',[num2str(Angle(1)),' ',num2str(Angle(2))]);
set(sdip.hdl.int,'String',Ijs_m);

%________________________________________________________________________
case 'clearall'
%------------------------------------------------------------------------
% Clears all dipoles, and reset the toggle buttons

if isfield(st.vols{1},'sdip')
    sdip = st.vols{1}.sdip;
else
    error('I can''t find sdip structure');
end

disp('Clears all dipoles')
ImaGIN_DrawElectrodes('ClearDip');
for ii=1:st.vols{1}.sdip.n_seeds
    if sdip.exitflag(ii)==1
        set(st.vols{1}.sdip.hdl.hseed(ii),'Value',0);
    end
end
set(st.vols{1}.sdip.hdl.hdip,'Value',1);


%________________________________________________________________________
case 'chgdip'
%------------------------------------------------------------------------
% Changes the dipole index for the first seed displayed

sdip = st.vols{1}.sdip;

i_dip = get(sdip.hdl.hdip,'Value');
if gco==sdip.hdl.hdipplus
    i_dip=i_dip+1;
    set(sdip.hdl.hdip,'Value',i_dip)
elseif gco==sdip.hdl.hdipminus
    i_dip=i_dip-1;
    set(sdip.hdl.hdip,'Value',i_dip)
end
if i_dip>=sdip.n_dip
    set(sdip.hdl.hdipplus,'Visible','off')
else
    set(sdip.hdl.hdipplus,'Visible','on')    
end
if i_dip<=1
    set(sdip.hdl.hdipminus,'Visible','off')
else
    set(sdip.hdl.hdipminus,'Visible','on')    
end
i_seed=1;
if ~isempty(i_seed)
    ImaGIN_DrawElectrodes('ClearDip')
    ImaGIN_DrawElectrodes('DrawDip',i_seed,i_dip);
end

%________________________________________________________________________
case 'cleardip'
%------------------------------------------------------------------------
% FORMAT ImaGIN_DrawElectrodes('ClearDip',seed_i)
% e.g. ImaGIN_DrawElectrodes('ClearDip')
%       clears all displayed dipoles
% e.g. ImaGIN_DrawElectrodes('ClearDip',1)
%       clears the first dipole displayed

if nargin>2
    seed_i = varargin{1};
else
    seed_i = 0;
end
if isfield(st.vols{1},'sdip')
        sdip = st.vols{1}.sdip;
    else
        return; % I don't do anything, as I can't find sdip strucure
end

if isfield(sdip,'ax')
    Nax = size(sdip.ax,2);
    else
        return; % I don't do anything, as I can't find axes info
end

if seed_i==0 % removes everything
    for ii=1:Nax
        for jj=1:6
            delete(sdip.ax(jj,ii));
        end
    end
    sdip = rmfield(sdip,'tabl_seed_dip');
    sdip = rmfield(sdip,'ax');    
elseif seed_i<=Nax % remove one seed only
    l_seed = find(sdip.tabl_seed_dip(:,1)==seed_i);
    for ii=l_seed
        for jj=1:6
            delete(sdip.ax(jj,ii));
        end
    end
    sdip.ax(:,l_seed) = [];
    sdip.tabl_seed_dip(l_seed,:) = [];
else
    error('Trying to clear unspecified dipole');
end
st.vols{1}.sdip = sdip;

%________________________________________________________________________
otherwise,
	warning('Unknown action string')
end;

end

%________________________________________________________________________
%________________________________________________________________________
%
% SUBFUNCTIONS
%________________________________________________________________________
%________________________________________________________________________
function dh = add1dip(loc,js,mark,col,ax,Fig,bb)
% Plots the dipoles on the 3 views
% Then returns the handle to the plots

loc(1,:) = loc(1,:) - bb(1,1)+1;
loc(2,:) = loc(2,:) - bb(1,2)+1;
loc(3,:) = loc(3,:) - bb(1,3)+1;
% +1 added to be like John's orthview code

dh = zeros(6,1);
figure(Fig)
% Transverse slice, # 1
set(Fig,'CurrentAxes',ax{1}.ax)
set(ax{1}.ax,'NextPlot','add')
dh(1) = plot(loc(1),loc(2),[mark,col],'LineWidth',2);
dh(2) = plot(loc(1)+[0 js(1)],loc(2)+[0 js(2)],col,'LineWidth',2);
set(ax{1}.ax,'NextPlot','replace')

% Coronal slice, # 2
set(Fig,'CurrentAxes',ax{2}.ax)
set(ax{2}.ax,'NextPlot','add')
dh(3) = plot(loc(1),loc(3),[mark,col],'LineWidth',2);
dh(4) = plot(loc(1)+[0 js(1)],loc(3)+[0 js(3)],col,'LineWidth',2);
set(ax{2}.ax,'NextPlot','replace')

% Sagital slice, # 3
set(Fig,'CurrentAxes',ax{3}.ax)
set(ax{3}.ax,'NextPlot','add')
dh(5) = plot(bb(2,2)-bb(1,2)-loc(2),loc(3),[mark,col],'LineWidth',2);
dh(6) = plot(bb(2,2)-bb(1,2)-loc(2)+[0 -js(2)],loc(3)+[0 js(3)],col,'LineWidth',2);
set(ax{3}.ax,'NextPlot','replace')

end


