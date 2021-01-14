function [nframes,res,fov,timeres,VENC,area_vol,diam_vol,flowPerHeartCycle_vol, ...
    maxVel_vol,PI_vol,RI_vol,flowPulsatile_vol,velMean_val, ...
    VplanesAllx,VplanesAlly,VplanesAllz,Planes,branchList,segment,r, ...
    timeMIPcrossection,segment1,vTimeFrameave,MAGcrossection, ...
    bnumMeanFlow,bnumStdvFlow,StdvFromMean] = loadpcvipr(directory,handles)
%LOADPCVIPR: loadpcvipr reads in header information and reconstructed data 
%(velocity, vmean, etc.) and transforms data into usable matlab variables.
%
%   Used by: paramMap.m
%   Dependencies: load_dat.m, background_phase_correction.m, evaluate_poly.m
%     calc_angio.m, feature_extraction.m, paramMap_params_new.m, makeITPlane.m
%     slidingThreshold.m

%% Reads PCVIPR Header
fid = fopen([directory '\pcvipr_header.txt'], 'r');
delimiter = ' ';
formatSpec = '%s%s%[^\n\r]'; %read 2 strings(%s%s),end line(^\n),new row(r)
% Info from headers are placed in dataArray, 1x2 cell array.
dataArray = textscan(fid, formatSpec, 'Delimiter', delimiter, ...
    'MultipleDelimsAsOne', true, 'ReturnOnError', false);
fclose(fid);

% Converts value column from strings to structure with nums.
dataArray{1,2} = cellfun(@str2num,dataArray{1,2}(:), 'UniformOutput', false);
pcviprHeader = cell2struct(dataArray{1,2}(:), dataArray{1,1}(:), 1);

nframes = pcviprHeader.frames; %number of reconstructed frames
timeres = pcviprHeader.timeres; %temporal resolution
fov = (pcviprHeader.fovx)/10;  %field of view (cm)
res = pcviprHeader.matrixx; %number of pixels in rows (ASSUMED ISOTROPIC)
VENC = pcviprHeader.VENC;

if ~(res == pcviprHeader.matrixy && res == pcviprHeader.matrixz) 
    disp('Resolution is not isotropic!')
end 

%% Reads Data Header
% Checks if automatic background phase correction was performed in recon
fid = fopen([directory '\data_header.txt'], 'r');
if fid>0
    dataHeader = textscan(fid, formatSpec, 'Delimiter', delimiter, ...
        'MultipleDelimsAsOne', true, 'ReturnOnError', false);
    fclose(fid);
    bgpcIdx = find(contains(dataHeader{1,1},'automatic_BGPC_flag'));
    if isempty(bgpcIdx)
        BGPCdone = 0;
    else
        BGPCdone = dataHeader{1,2}{bgpcIdx};
        BGPCdone = str2double(BGPCdone);
    end 
else
    BGPCdone = 0; %assume automatic backg. phase corr. wasnt done in recon
end 

%% Read MAG Data    
set(handles.TextUpdate,'String','Loading Time Averaged Data'); drawnow;
MAG = load_dat(fullfile(directory,'MAG.dat'),[res res res]);
  
%% Auto crop images (from MAG data)
% Done to save memory when loading in TR velocity data below.
SUMnumA = squeeze(sum(sum(MAG,1),2)); %1D axial projection
SUMnumS = squeeze(sum(sum(MAG,1),3))'; %1D sagittal projection
SUMnumC = squeeze(sum(sum(MAG,2),3)); %1D coronal projection

% Chop off edges of projections (usually noisy data)
SUMnumA(1:3) = 0; 
SUMnumA(end-2:end) = 0;
SUMnumS(1:3) = 0;
SUMnumS(end-2:end) = 0;
SUMnumC(1:3) = 0;
SUMnumC(end-2:end) = 0;

% Combine projections and normalize values (from 0-1)
SUMnum = [SUMnumC,SUMnumS,SUMnumA];
SUMnum = rescale(SUMnum,'InputMin',min(SUMnum),'InputMax',max(SUMnum));

% Find where projection crosses normalized threshold value of 0.25 
BIN = SUMnum>0.25;
[~,IDXstart] = max(BIN,[],1); %get first thresh crossing
[~,IDXend] = max(flipud(BIN),[],1); 
IDXend = size(MAG) - IDXend + 1; %get last thresh crossing

% Crop data with new dimensions
MAG = MAG(IDXstart(1):IDXend(1),IDXstart(2):IDXend(2),IDXstart(3):IDXend(3));
newDIM = size(MAG);
    
%% Read Average Velocity
vMean = zeros(newDIM(1),newDIM(2),newDIM(3),3,'single'); 
% Looped reading of average velocity data.
for n = 1:3
    temp = load_dat(fullfile(directory,['comp_vd_' num2str(n) '.dat']),[res res res]);
    vMean(:,:,:,n) = temp(IDXstart(1):IDXend(1),IDXstart(2):IDXend(2),IDXstart(3):IDXend(3));
end

%% Create Angio
% Calculate complex difference angiogram for visualization.
set(handles.TextUpdate,'String','Creating Angiogram'); drawnow;
timeMIP = calc_angio(MAG, vMean, VENC);
% NOTE: timeMIP is an approximated complex difference image.
% The result is nearly equivalent to loading 'CD.dat'.

%% Find optimum global threshold 
step = 0.001; %step size for sliding threshold
UPthresh = 0.8; %max upper threshold when creating Sval curvature plot
SMf = 10;
shiftHM_flag = 1; %flag to shift max curvature by FWHM
[~,segment] = slidingThreshold(timeMIP,step,UPthresh,SMf,shiftHM_flag);
areaThresh = round(sum(segment(:)).*0.005); %minimum area to keep
conn = 6; %connectivity (i.e. 6-pt)
segment = bwareaopen(segment,areaThresh,conn); %inverse fill holes

clear CDcrop x y SMf temp n halfMaxRightIndex halfMaxLeftIndex Idx BIN 
clear curvatureSM denom num ddy dy ddx dx areaThresh fullWidth conn
clear Sval iter maxThresh newDIM dataArray fid formatSpec delimiter
clear SUMnumA SUMnumC SUMnumS SUMnum step ans dataHeader UPthresh shiftHM_flag

%% Feature Extraction
% Get trim and create the centerline data
sortingCriteria = 3; %sorts branches by junctions/intersects 
spurLength = 15; %minimum branch length (removes short spurs)
[~,~,branchList,~] = feature_extraction(sortingCriteria,spurLength,vMean,segment,handles);

% Flow parameter calculation, bulk of code is in paramMap_parameters.m
[area_vol,diam_vol,flowPerHeartCycle_vol,maxVel_vol,PI_vol,RI_vol, ...
    flowPulsatile_vol,velMean_val,VplanesAllx,VplanesAlly,VplanesAllz, ... 
    r,timeMIPcrossection,segment1,vTimeFrameave,MAGcrossection, ...
    bnumMeanFlow,bnumStdvFlow,StdvFromMean,Planes] ...
    = paramMap_params_new(branchList,res,timeMIP,vMean,directory, ...
    BGPCdone,nframes,fov,MAG,IDXstart,IDXend,handles);

set(handles.TextUpdate,'String','All Data Loaded'); drawnow;

return