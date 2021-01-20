% This script reads and displays gating track information.

% Author: Grant Roberts
% Date: September 18 2019


%% Filter Gating Directory
gatingDir = uigetdir(pwd, 'Select the directory with the gating tracks');
cd(gatingDir) % move to this directory

fullDir = dir; % get all files in this directory
% remove all checksum and h5 files
keepFiles = ~(contains({fullDir.name},'.md5sum') | contains({fullDir.name},'.h5'))';
filtDir = fullDir(keepFiles); % filter out md5sum and h5 files
fileNames = {filtDir.name}'; % get filenames from filtered directory
filesLoaded = {}; % create a cell to add loaded file names
ecgData.gatingDir = gatingDir; % add gating directory to a composite data structure
diary ecgInformation.txt % start writing command window output to text file



%% Gating Tracks
badProspGating=0;

%Load Full Gating Track - Shows full gating information
idx = find(contains(fileNames,'pcvipr_track.full')); % does this file exist in fileNames?
if ~isempty(idx) % if we found something..
    name = fileNames{idx}; % grab name from fileNames
    fid = fopen(name);
    gate = fread(fid,'int32','b');
    gate = reshape(gate,[numel(gate)/5 5]); % data is Nx5 format
    gate = sortrows(gate,3); % organize time chronologically (if not already)
    fclose(fid);

    % Put data into gatingTrackFull structure
    gatingTrackFull.ecg = gate(:,1); % ecg (sawtooth) waveform
    gatingTrackFull.resp = 4095-gate(:,2); % resp waveform (12bit max)
    gatingTrackFull.time = gate(:,3)/1e3; % time at each data point (ms)
    gatingTrackFull.prep = gate(:,4); % unsure what this is
    gatingTrackFull.acq = gate(:,5); % data point acquisition number
    [rrFull, peaksFull] = getRR(gatingTrackFull.ecg); % get all rr intervals
        gatingTrackFull.rr = (nan(size(gate,1),1)); % expand rr vector to match gate lengths
        gatingTrackFull.rr(peaksFull) = rrFull; % fill in expanded Nan vector with rrs
    encodeLength = mode(hist(gatingTrackFull.acq,unique(gatingTrackFull.acq))); % get encoding type (5-point e.g.)
        gatingTrack.encodeLength = encodeLength;
    [missedHBidx,earlyTrigIdx] = findMissedHB(rrFull,peaksFull); % detect missed heartbeats and early triggers
        numMissedFull = length(missedHBidx);
        numEarlyFull = length(earlyTrigIdx);
        missedHBs = nan(length(gatingTrackFull.rr),1); % expand vector to match gate lengths
        earlyTrigs = nan(length(gatingTrackFull.rr),1); % expand vector to match gate lengths
        missedHBs(missedHBidx) = gatingTrackFull.rr(missedHBidx); % fill in expanded vector with missed HBs
        earlyTrigs(earlyTrigIdx) = gatingTrackFull.rr(earlyTrigIdx); % fill in expanded vector early triggers
        gatingTrackFull.missed = missedHBs;
        gatingTrackFull.early = earlyTrigs;
        gatingTrackFull.clean = gatingTrackFull.rr;
        gatingTrackFull.clean(missedHBidx) = nan;
        gatingTrackFull.clean(earlyTrigIdx) = nan; % data which has missed HB and early trigger events removed
        gatingTrackFull.bpm = 60000./(gatingTrackFull.clean); % convert rr (ms) to heartrate (bpm)
              
    fileNames(idx) = []; % clear name from filName list to allow next gating track to be loaded.
    filesLoaded = [filesLoaded, name]; % add file to loaded files running cell
    ecgData.gatingTrackFull = gatingTrackFull; % composite structure containing all information
end

%Load Gating Track - Shows nly prospectively acquired gating data
idx = find(contains(fileNames,'pcvipr_track'));
if ~isempty(idx)      
    name = fileNames{idx};
    fid = fopen(name);
    gate = fread(fid,'int32','b');
    gate = reshape(gate,[numel(gate)/4 4]);
    gate = sortrows(gate,3);
    fclose(fid);

    % Put data into gatingTrack structure
    gatingTrack.ecg = gate(:,1);
    gatingTrack.resp = 4095-gate(:,2);
    gatingTrack.time = gate(:,3)/1e3;
    gatingTrack.prep = gate(:,4);
    gatingTrack.acq = (1:size(gate,1))'; % generate a data acquisition row
    [rr, peaks] = getRR(gatingTrack.ecg);
        gatingTrack.rr = (nan(size(gate,1),1));
        gatingTrack.rr(peaks) = rr;
    [missedHBidx,earlyTrigIdx] = findMissedHB(rr,peaks);
        numMissed = length(missedHBidx);
        numEarly = length(earlyTrigIdx);
        missedHBs = nan(length(gatingTrack.rr),1);
        earlyTrigs = nan(length(gatingTrack.rr),1);
        missedHBs(missedHBidx) = gatingTrack.rr(missedHBidx);
        earlyTrigs(earlyTrigIdx) = gatingTrack.rr(earlyTrigIdx);
        gatingTrack.missed = missedHBs;
        gatingTrack.early = earlyTrigs;
        gatingTrack.clean = gatingTrack.rr;
        gatingTrack.clean(missedHBidx) = nan;
        gatingTrack.clean(earlyTrigIdx) = nan;
        gatingTrack.bpm = 60000./gatingTrack.clean;
    filesLoaded = [filesLoaded;name];
    ecgData.gatingTrack = gatingTrack;
    if mode(gatingTrack.time)==0 % if the data is not being acquired for most of the scan
        fprintf('Error during acqusition.\n'); % then we have an issue with the prospective gating.
        fprintf('The pcvipr_track file shows that data was not being acquired consistently. Script may crash.\n');
        return % kick us out
    end 
end

if sum(contains(filesLoaded,'pcvipr_track'))==2 % if we have both ".pcvipr_track" and ".pcvipr_track.full" files..
    if gatingTrackFull.time(end) < gatingTrack.time(end) % if the end times don't match up..
        badProspGating=1; % then the prospective gating was bad.
        fprintf('The .pcvipr_track file shows that the time of acquisition is larger than the time of the full waveform.\n')
        fprintf('This is likely an acquisition an error from dysfunctional prospective gating. Script may crash.\n\n');
    end 
end 



%% ECG Track (old gating)
ecgIdx = find(contains(fileNames,'ecg_track')); % does this file exist in fileNames?
if ~isempty(ecgIdx) % if we found something..
    ecgName = fileNames{ecgIdx}; % grab name from fileNames
    infoIdx = find(contains(fileNames,'scan_info.txt')); % do we have the scan_info.txt file?
    if isempty(infoIdx) % if we found scan_info.txt ..
        fprintf('  \t A "scan_info.txt" file is required to load this ecg_track file\n');
        exit
    else 
        infoName = fileNames{infoIdx}; % grab name from fileNames
        fid = fopen(infoName);  % open scan_info.txt

        needTR = 1; % flag if we still need to find TR
        needNPROJ = 1; % flag if we still need to find # projections
        needINTER = 1; % flag if we still need to find interleave #
        needORDER = 1; % flag if we still need to find projection order
        looking = needTR + needNPROJ + needINTER + needORDER; % if > 0, keep looking
        count = 0;
        while looking                                     
            tline = fgetl(fid); % grab individual lines from scan_info (no easy way to do this)                             
            if ischar(tline)                       
                if needTR
                    param1 = contains(tline, 'TR'); % find the "TR" line  
                    if param1             
                        needTR = 0; % turn flag off to stop searching
                        %tr=sscanf(tline, '%*s %*s %s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s');
                        tr = sscanf(tline, '%*s %*s %f');
                        %tr=str2double(tr); % turn from string into number
                    end 
                end 
                if needNPROJ
                    param2 = contains(tline,'rhuser8') | contains(tline,'rhuser08'); 
                    if param2                
                        needNPROJ = 0; 
                        %nproj=sscanf(tline, '%*s %*s %s');
                        nproj=sscanf(tline, '%*s %*s %d');
                        %nproj=str2double(nproj);
                    end 
                end 
                if needINTER
                    param3 = contains(tline, 'rhuser10'); 
                    if param3            
                        needINTER = 0; 
                        inter = sscanf(tline, '%*s %*s %d');
                    end 
                end 
                if needORDER
                    param4 = strfind(tline, 'PR Order');
                    if isfinite(param4) == 1                
                        needORDER = 0; 
                        PR_ORDER = sscanf(tline, ' %*s %*s %s');
                        PR_ORDER = str2double(PR_ORDER);
                    end 
                end   
            end
            looking = needTR + needNPROJ + needINTER + needORDER; % update while loop flag
            if count>3000
                break
            end 
            count = count+1;
            PR_ORDER = 2;
        end  

        % Load ECG Track
        subproj = nproj / inter; % calculate # projections in interleave
        fid = fopen(ecgName, 'rb');
        ecg = fread(fid,'int', 'b');
        fclose(fid);

        % Determine projection order, changes depending on PR_ORDER
        acquisition_order = zeros(nproj,1); % array for all projections     
        sub_order = bitreverse(subproj); % bit reverse for the subframes

        if PR_ORDER == 2 
            frame_order = 0:inter-1; % subframes are ordered sequentially
        elseif PR_ORDER == 4
            frame_order = bitrevorder(inter); % subframes are bit reversed
        end

        for pos = 0:nproj-1   
            inter_index = mod(pos,inter); % which interleave am I in?
            inter_pos = floor(pos/inter); % index position w/in interleave
            sub_pos = inter_index + inter_pos; % not sure what this is
            if PR_ORDER == 2
                acquisition_order(pos+1) = subproj*frame_order(floor(pos/subproj)+1) + sub_order(mod(pos,subproj)+1);
            end

            if  PR_ORDER == 4
                acquisition_order(pos+1) = sub_order(mod(sub_pos,subproj)+1) + subproj*frame_order(inter_index+1);
            end
        end

        ecg_sorted = ecg(acquisition_order+1); %amplitude
        TR=tr*4; % for 4point encode
        time = (TR:TR:TR*length(ecg_sorted))'; %time in ms
    end 

    % Put data into ecgTrack structure
    ecgTrack.ecg = ecg_sorted; % ecg (sawtooth) waveform
    ecgTrack.time = time; % time at each data point (ms)
    ecgTrack.acq = (1:length(ecg_sorted))'; % data acquisition number
    [rrECG, peaksECG] = getRR(ecgTrack.ecg); % get rr intervals from ecg data
        ecgTrack.rr = nan(length(ecgTrack.ecg),1); % expand rr vector to match length of ecg vector
        ecgTrack.rr(peaksECG) = rrECG; % fill in NaNs with rr intervals
    [missedHBidx,earlyTrigIdx] = findMissedHB(rrECG,peaksECG); % find missed HBs and early triggers
        numMissedECG = length(missedHBidx);
        numEarlyECG = length(earlyTrigIdx);
        missedHBs = nan(length(ecgTrack.rr),1); % expand vector to match ecg length
        earlyTrigs = nan(length(ecgTrack.rr),1); % expand vector to match ecg length
        missedHBs(missedHBidx) = ecgTrack.rr(missedHBidx); % fill in vector with bad data from rr 
        earlyTrigs(earlyTrigIdx) = ecgTrack.rr(earlyTrigIdx); % fill in vector with bad data from rr
        ecgTrack.missed = missedHBs;
        ecgTrack.early = earlyTrigs;
        ecgTrack.clean = ecgTrack.rr;
        ecgTrack.clean(missedHBidx) = nan;
        ecgTrack.clean(earlyTrigIdx) = nan; % rr data without missed HBs and early triggers
        ecgTrack.bpm = 60000./(ecgTrack.clean); % convert rr (ms) to heartrate (bpm)
    filesLoaded = [filesLoaded;ecgName]; % add file to loaded files running cell

end



%% Complete ECG, PPG, Respiratory Waveforms

% Second ECG Lead - Full Waveform. If PPG gated, this will be noise.
idx = find(contains(fileNames,'ECG2Data')); % does this file exist in fileNames?
if ~isempty(idx) % if we found something..
    name = fileNames{idx}; % grab name from fileNames
    waveform.ecg2 = importdata(name); % put data into waveform structure
    filesLoaded = [filesLoaded;name];
    ecgData.waveform.ecg2 = waveform.ecg2;
end 

% Third ECG Lead - Full Waveform. If PPG gated, this will be noise.
idx = find(contains(fileNames,'ECG3Data'));
if ~isempty(idx)
    name = fileNames{idx};
    waveform.ecg3 = importdata(name);
    filesLoaded = [filesLoaded;name];
    ecgData.waveform.ecg3 = waveform.ecg3;
end 

% Peripheral Gating - Full Waveform. % If ECG gated, this will be noise.
% Obtained via pulse oximetry on index finger.
idx = find(contains(fileNames,'PPGData'));
if ~isempty(idx)
    name = fileNames{idx};
    waveform.ppg = importdata(name);
    filesLoaded = [filesLoaded;name];
    ecgData.waveform.ppg = waveform.ppg;
end 

% Respiratory Gating - Full Waveform. Taken from respiratory bellows.
idx = find(contains(fileNames,'RESPData'));
if ~isempty(idx)
    name = fileNames{idx};
    waveform.resp = importdata(name);
    waveform.resp = 4095-waveform.resp;
    filesLoaded = [filesLoaded;name];
    ecgData.waveform.resp = waveform.resp;
end 



%% Recorded ECG, PPG, Respiratory Triggers

% QRS Triggers from Second ECG Lead. If PPG gated, this will be noise.
idx = find(contains(fileNames,'ECG2Trig'));  % does this file exist in fileNames?
if ~isempty(idx) % if we found something..
    name = fileNames{idx}; % grab name from fileNames
    trigger.ecg2 = importdata(name); % put data into trigger structure
    filesLoaded = [filesLoaded;name];
    ecgData.trigger.ecg2 = trigger.ecg2;
end 

% QRS Triggers from Third ECG Lead. If PPG gated, this will be noise.
idx = find(contains(fileNames,'ECG3Trig'));
if ~isempty(idx)
    name = fileNames{idx};
    trigger.ecg3 = importdata(name);
    filesLoaded = [filesLoaded;name];
    ecgData.trigger.ecg3 = trigger.ecg3;
end 

% Peak Signal Trigger from Peripheral Gating.
% Obtained via pulse oximetry on index finger.
idx = find(contains(fileNames,'PPGTrig'));
if ~isempty(idx)
    name = fileNames{idx};
    trigger.ppg = importdata(name);
    filesLoaded = [filesLoaded;name];
    ecgData.trigger.ppg = trigger.ppg;
end 

% Triggers from Respiratory Bellows
idx = find(contains(fileNames,'RESPTrig'));
if ~isempty(idx)
    name = fileNames{idx};
    trigger.resp = importdata(name);
    filesLoaded = [filesLoaded;name];
    ecgData.trigger.resp = trigger.resp;
end 



%% Plotting (just gating track, no full)
if sum(contains(filesLoaded,'pcvipr_track'))==1 || badProspGating % only the "pcvipr_track" file was loaded or we have bad gating..
    dt = diff(gatingTrack.time); % take derivative of time vector
    timeres = mode(dt); % find time resolution by finding most common value
    iterator = round(gatingTrack.time/timeres); % get number of "should be" data points. gatingTrack doesn't include data if acquisition isn't on.
    respOn = nan(iterator(end),1); % create expanded resp vector
    gateOn = nan(iterator(end),1); % create expanded rr interval vector
    expandedTime = (0:timeres:gatingTrack.time(end))'; % create expanded time vector, filling in missing timepoints

    for i=1:length(iterator)
        respOn(iterator(i)) = gatingTrack.resp(i); % place known respiratory data in expanded respOn vector
        gateOn(iterator(i)) = gatingTrack.rr(i); % place known ecg data in expanded gateOn vector
    end 
    gatingTrack.timeres = timeres; % save time resolution
    gatingTrack.respOn = respOn; % save respiratory on time
    gatingTrack.gateOn = gateOn; % save ecg on time

    % Respiratory Subfigure
    figure; hFig = figure(1); set(hFig, 'Position', [50 50 1800 800]); % enlarge figure
    sgtitle('Gating Track Information'); % create overall title
    subplot(2,1,1); hold on; title('Respiratory Waveforms'); % make subplot and subplot title
    plot(expandedTime/1000,respOn,'Color',[0, 0.5, 0.19],'LineWidth',2); % plot expanded respiratory waveform
    xlabel('Time (s)'); ylabel('Respiratory Amplitude (a.u.)');
    legend('Respiratory Gating');

    % PG/ECG Subfigure
    maxRR = ones(length(expandedTime/1000),1).*nanmax(gatingTrack.clean); % get maximum rr interval from clean RRs
    meanRR = ones(length(expandedTime/1000),1).*nanmean(gatingTrack.clean); % get mean rr interval
    minRR = ones(length(expandedTime/1000),1).*nanmin(gatingTrack.clean); % get mininimum rr interval

    subplot(2,1,2); hold on; title('ECG/PG Triggers (BPM)') % make second subplot
    plot(expandedTime/1000,maxRR,'--','Color','black'); % plot dashed line of max RR
    plot(expandedTime/1000,meanRR,'--','Color','black','HandleVisibility','off'); % plot dashed line of mean RR
    plot(expandedTime/1000,minRR,'--','Color','black','HandleVisibility','off'); % plot dashed line of min RR
    scatter(expandedTime/1000,gateOn,35,'blue'); % scatter plot of rr intervals over time
    scatter(gatingTrack.time/1000,gatingTrack.missed,35,'filled','red'); % fill in points with red if missed HB
    scatter(gatingTrack.time/1000,gatingTrack.early,35,'filled','yellow'); % fill in points with yellow if early trig
    xlabel('Time (s)'); ylabel('RR Interval (ms)');
    legend('Max/Min/Mean RR','ECG/PG RR-intervals','Likely Missed HBs','Likely Early Triggers');
    savefig('gatingTrack_resp_ecg_waveforms'); % save this figure
    
    % Histogram (BPM)
    figure; histogram(gatingTrack.rr); % plot histogram of rr intervals (not clean rr intervals)
    xlabel('RR-interval (ms)'); ylabel('Frequency (counts)'); title('Histogram of RR intervals');
    savefig('gatingTrack_histogram'); % save this figure
    
    
    %%%%% Command Line Output %%%%%
    fprintf('ECG/PG Gating Information and Statistics:\n') 
    fprintf('\n');
    fprintf('    Time: \t \t \t \t %s\n',datetime('now'));
    fprintf('    Directory: \t \t \t %s\n',gatingDir);
    for i=1:length(filesLoaded)
        fprintf('    Files Loaded: \t \t %s\n',filesLoaded{i});
    end 
    fprintf('\n');
    if length(expandedTime)>length(gatingTrack.time)
        fprintf('    Gating Type: \t \t \t \t \t \t Prospective Gating\n');
        fprintf('    Percent Data Acquired: \t \t \t \t %7.2f %%\n', 100*length(gatingTrack.time)/length(expandedTime));    
    else
        fprintf('    Gating Type: \t \t \t \t \t Retrospective Gating\n');
    end 
    fprintf('    ECG/Resp. time resolution: \t \t \t %7.2f ms.\n',gatingTrack.timeres);
    fprintf('    Duration: \t \t \t \t \t \t \t %7.2f s.\n',gatingTrack.time(end)/1000);
    fprintf('\n');
    fprintf('    Mean RR interval: \t \t \t \t \t %7.2f ms.\n',nanmean(gatingTrack.clean));
    fprintf('    Median RR interval: \t \t \t \t %7.2f ms.\n',nanmedian(gatingTrack.clean));
    fprintf('    Standard deviation: \t \t \t \t %7.2f ms.\n',nanstd(gatingTrack.clean));
    fprintf('    Range: \t \t \t \t \t \t \t \t %7.2f ms.\n',range(gatingTrack.clean));
    fprintf('\n');
    fprintf('    Mean heart rate: \t \t \t \t \t %7.2f bpm.\n',nanmean(gatingTrack.bpm));
    fprintf('    Median heart rate: \t \t \t \t \t %7.2f bpm.\n',nanmedian(gatingTrack.bpm));
    fprintf('    Standard deviation: \t \t \t \t %7.2f bpm.\n',nanstd(gatingTrack.bpm));
    fprintf('    Range: \t \t \t \t \t \t \t \t %7.2f bpm.\n',range(gatingTrack.bpm));
    fprintf('    Total Recorded heartbeats: \t \t \t \t %d beats.\n',length(peaks));
    fprintf('    Estimated # missed heartbeats: \t \t \t \t %d beats (%2.2f%%).\n',numMissed,100*(numMissed/length(peaks)));
    fprintf('    Estimated # early triggers: \t \t \t \t %d beats (%2.2f%%).\n\n',numEarly,100*(numEarly/length(peaks)));
end 
        
        
        
%% Plotting (ecg track)
if contains(filesLoaded,'ecg_track') % if we loaded in an ecg_track file ..
    dt = diff(ecgTrack.time);
    timeres = mode(dt);
    iterator = round(ecgTrack.time/timeres);
    gateOn = nan(iterator(end),1);
    expandedTime = (0:timeres:ecgTrack.time(end))';

    for i=1:length(iterator)
        gateOn(iterator(i)) = ecgTrack.rr(i);
    end 
    ecgTrack.timeres = timeres;
    ecgTrack.gateOn = gateOn;
    
    % PG figure
    maxRR = ones(length(expandedTime/1000),1).*nanmax(ecgTrack.clean);
    meanRR = ones(length(expandedTime/1000),1).*nanmean(ecgTrack.clean);
    minRR = ones(length(expandedTime/1000),1).*nanmin(ecgTrack.clean);
    
    figure; hFig = figure(1); set(hFig, 'Position', [50 50 1800 800]); 
    hold on; title('ECG/PG Triggers (BPM)')
    plot(expandedTime/1000,maxRR,'--','Color','black');
    plot(expandedTime/1000,meanRR,'--','Color','black','HandleVisibility','off');
    plot(expandedTime/1000,minRR,'--','Color','black','HandleVisibility','off');
    scatter(expandedTime/1000,gateOn,35,'blue');
    scatter(ecgTrack.time/1000,ecgTrack.missed,35,'filled','red');
    scatter(ecgTrack.time/1000,ecgTrack.early,35,'filled','yellow');
    xlabel('Time (s)'); ylabel('RR Interval (ms)');
    legend('Max/Min/Mean RR','ECG/PG RR-intervals','Likely Missed HBs','Likely Early Triggers');
    savefig('ecgTrack_waveform');
    
    % Histogram (BPM)
    figure; histogram(ecgTrack.rr);
    xlabel('RR-interval (ms)'); ylabel('Frequency (counts)'); title('Histogram of RR intervals');
    savefig('ecgTrack_histogram');
    
    %%%%% Command Line Output %%%%%
    fprintf('ECG/PG Gating Information and Statistics:\n')
    fprintf('\n');
    fprintf('    Time: \t \t \t \t %s\n',datetime('now'));
    fprintf('    Directory: \t \t \t %s\n',gatingDir);
    for i=1:length(filesLoaded)
        fprintf('    Files Loaded: \t \t %s\n',filesLoaded{i});
    end 
    fprintf('\n');
    if length(expandedTime)>length(ecgTrack.time)
        fprintf('    Gating Type: \t \t \t \t \t \t Prospective Gating\n');
        fprintf('    Percent Data Acquired: \t \t \t \t %7.2f %%\n', 100*length(ecgTrack.time)/length(expandedTime));    
    else
        fprintf('    Gating Type: \t \t \t \t \t Retrospective Gating\n');
    end 
    fprintf('    ECG time resolution: \t \t %7.2f ms.\n',ecgTrack.timeres);
    fprintf('    Duration: \t \t \t \t \t \t \t %7.2f s.\n',ecgTrack.time(end)/1000);
    fprintf('\n');
    fprintf('    Mean RR interval: \t \t \t \t \t %7.2f ms.\n',nanmean(ecgTrack.clean));
    fprintf('    Median RR interval: \t \t \t \t %7.2f ms.\n',nanmedian(ecgTrack.clean));
    fprintf('    Standard deviation: \t \t \t \t %7.2f ms.\n',nanstd(ecgTrack.clean));
    fprintf('    Range: \t \t \t \t \t \t \t \t %7.2f ms.\n',range(ecgTrack.clean));
    fprintf('\n');
    fprintf('    Mean heart rate: \t \t \t \t \t %7.2f bpm.\n',nanmean(ecgTrack.bpm));
    fprintf('    Median heart rate: \t \t \t \t \t %7.2f bpm.\n',nanmedian(ecgTrack.bpm));
    fprintf('    Standard deviation: \t \t \t \t %7.2f bpm.\n',nanstd(ecgTrack.bpm));
    fprintf('    Range: \t \t \t \t \t \t \t \t %7.2f bpm.\n',range(ecgTrack.bpm));
    fprintf('    Estimated # missed heartbeats: \t \t \t \t %d beats (%2.2f%%).\n',numMissedECG,100*(numMissedECG/length(peaksECG)));
    fprintf('    Estimated # early triggers: \t \t \t \t %d beats (%2.2f%%).\n\n',numEarlyECG,100*(numEarlyECG/length(peaksECG))); 
end 



%% Plotting (with full gating track)
if sum(contains(filesLoaded,'pcvipr_track'))==2 % if we have both ".pcvipr_track" and ".pcvipr_track.full" files..
    if mod(length(gatingTrackFull.acq),length(gatingTrack.acq))==0 % if gatingTrackFull.acq is a multiple of the gatingTrack.acq
        prospectiveGating = 0; % set prospective gating flag to false
    else
        prospectiveGating = 1; % set flag to true
    end

    dt = diff(gatingTrack.time);
    timeres = mode(dt);
    timeresFull = mode(diff(gatingTrackFull.time));
    iterator = round(gatingTrack.time/timeres);
    respOn = nan(iterator(end),1);
    gateOn = nan(iterator(end),1);
    expandedTime = gatingTrackFull.time(encodeLength:encodeLength:end);

    for i=1:length(iterator)
        respOn(iterator(i)) = gatingTrack.resp(i);
        gateOn(iterator(i)) = gatingTrack.rr(i);
    end 

    gatingTrack.timeres = timeres;
    gatingTrackFull.timeres = timeresFull;
    gatingTrack.respOn = respOn;
    gatingTrack.gateOn = gateOn;

    if ~(length(respOn)==length(expandedTime)) % if the expanded vectors are not equal to the expanded time vectors..
        disp('Acquisition error. The pcvipr_track file shows that data was not being acquired consistently. Script may crash.');
        return % this happens if there are errors with scan acquisition
    end 


    % Respiratory Subfigure
    figure; hFig = figure(1); set(hFig, 'Position', [50 50 1800 800]);
    sgtitle('Gating Track and Full Gating Track Information'); 
    subplot(2,1,1); hold on; title('Respiratory Waveforms');
    plot(expandedTime/1000,respOn,'Color',[0, 0.5, 0.19],'LineWidth',2); 
    xlabel('Time (s)'); ylabel('Respiratory Amplitude (a.u.)');
    if prospectiveGating
        plot(gatingTrackFull.time/1000,gatingTrackFull.resp,'Color',[0.6, 0.6, 0.6]);
        legend('Respiratory Gating - Acquired Data','Respiratory Monitor (Full Waveform)');
    else 
        legend('Respiratory Gating');
    end 

    % PG Subfigure
    maxRR = ones(length(expandedTime/1000),1).*nanmax(gatingTrackFull.clean);
    meanRR = ones(length(expandedTime/1000),1).*nanmean(gatingTrackFull.clean);
    minRR = ones(length(expandedTime/1000),1).*nanmin(gatingTrackFull.clean);

    subplot(2,1,2); hold on; title('ECG/PG Triggers (BPM)')
    plot(expandedTime/1000,maxRR,'--','Color','black');
    plot(expandedTime/1000,meanRR,'--','Color','black','HandleVisibility','off');
    plot(expandedTime/1000,minRR,'--','Color','black','HandleVisibility','off');
    scatter(gatingTrackFull.time/1000,gatingTrackFull.rr,35,'blue');
    scatter(gatingTrackFull.time/1000,gatingTrackFull.missed,35,'filled','red');
    scatter(gatingTrackFull.time/1000,gatingTrackFull.early,35,'filled','yellow');
    xlabel('Time (s)'); ylabel('RR Interval (ms)');
    legend('Max/Min/Mean RR','ECG/PG RR-intervals','Likely Missed HBs','Likely Early Triggers');
    savefig('gatingTrack_resp_ecg_waveforms');

    % Histogram (BPM)
    figure; histogram(rrFull);
    xlabel('RR-interval (ms)'); ylabel('Frequency (counts)'); title('Histogram of RR intervals');
    savefig('gatingTrack_histogram');

    %%%%% Command Line Output %%%%%
    fprintf('ECG/PG Gating Information and Statistics:\n')
    fprintf('\n');
    fprintf('    Time: \t \t \t \t %s\n',datetime('now'));
    fprintf('    Directory: \t \t \t %s\n',gatingDir);
    for i=1:length(filesLoaded)
        fprintf('    Files Loaded: \t \t %s\n',filesLoaded{i});
    end 
    fprintf('\n');
    if length(expandedTime)>length(gatingTrack.time)
        fprintf('    Gating Type: \t \t \t \t \t \t Prospective Gating\n');
        fprintf('    Percent Data Acquired: \t \t \t \t %7.2f %%\n', 100*length(gatingTrack.time)/length(expandedTime));    
    else
        fprintf('    Gating Type: \t \t \t \t \t Retrospective Gating\n');
    end 
    fprintf('    ECG/Resp. time resolution (TR): \t %7.2f ms.\n',gatingTrackFull.timeres);
    fprintf('    Duration: \t \t \t \t \t \t \t %7.2f s.\n',gatingTrackFull.time(end)/1000);
    fprintf('\n');
    fprintf('    Mean RR interval: \t \t \t \t \t %7.2f ms.\n',nanmean(gatingTrackFull.clean));
    fprintf('    Median RR interval: \t \t \t \t %7.2f ms.\n',nanmedian(gatingTrackFull.clean));
    fprintf('    Standard deviation: \t \t \t \t %7.2f ms.\n',nanstd(gatingTrackFull.clean));
    fprintf('    Range: \t \t \t \t \t \t \t \t %7.2f ms.\n',range(gatingTrackFull.clean));
    fprintf('\n');
    fprintf('    Mean heart rate: \t \t \t \t \t %7.2f bpm.\n',nanmean(gatingTrackFull.bpm));
    fprintf('    Median heart rate: \t \t \t \t \t %7.2f bpm.\n',nanmedian(gatingTrackFull.bpm));
    fprintf('    Standard deviation: \t \t \t \t %7.2f bpm.\n',nanstd(gatingTrackFull.bpm));
    fprintf('    Range: \t \t \t \t \t \t \t \t %7.2f bpm.\n',range(gatingTrackFull.bpm));
    fprintf('    Total Recorded heartbeats: \t \t \t \t %d beats.\n',length(peaksFull));
    fprintf('    Estimated # missed heartbeats: \t \t \t \t %d beats (%2.2f%%).\n',numMissedFull,100*(numMissedFull/length(peaksFull)));
    fprintf('    Estimated # early triggers: \t \t \t \t %d beats (%2.2f%%).\n\n',numEarlyFull,100*(numEarlyFull/length(peaksFull)));
end 



%% Plotting (with complete waveforms)
if sum(contains(filesLoaded,'PPG'))==2 % if we have 2 ppg waveform files
    figure; plot(waveform.ppg); % unsure what the time resolution is..
    title('Raw PG Waveform with Triggers');
    hold on; trigPoints = nan(length(waveform.ppg),1); % expand trigger point vector to match ppg vector length
    trigPoints(trigger.ppg) = waveform.ppg(trigger.ppg); % find where the trigger points occur in waveform.ppg
    scatter(1:length(trigPoints),trigPoints,'filled','green'); % plot points of triggers on ppg waveform
    legend('Raw PPG Waveform','Triggers');
    fprintf('\nRaw Waveform Signal:\n\n');
    SNRppg = snr(waveform.ppg); % get SNR measure from power spectrum analysis
    if SNRppg>-8 % this is an arbitray threshold value, can be changed
        fprintf('   PPG: \t \t SIGNAL DETECTED\n');
    else
        fprintf('   PPG: \t \t LIKELY NOISE\n');
    end 
    savefig('Raw_PG_waveform_triggers.fig')
end 

if sum(contains(filesLoaded,'RESP'))==2 % if we have 2 resp waveform files
    figure; plot(waveform.resp);
    title('Raw Respiratory Waveform with Triggers');
    hold on; trigPoints = nan(length(waveform.resp),1);
    trigPoints(trigger.resp) = waveform.resp(trigger.resp);
    scatter(1:length(trigPoints),trigPoints,'filled','green');
    legend('Raw Respiratory Waveform','Bellow Triggers');
    SNRresp = snr(waveform.resp);
    if SNRresp>-8 % this is an arbitray threshold value, can be changed
        fprintf('   RESP: \t \t SIGNAL DETECTED\n');
    else
        fprintf('   RESP: \t \t LIKELY NOISE\n');
    end 
    savefig('Raw_resp_waveform_triggers.fig')
end 

if sum(contains(filesLoaded,'ECG'))==4 % if we have 4 ecg waveform files
    figure; sgtitle('Raw ECG Waveform with Triggers');
    subplot(2,1,1); plot(waveform.ecg2); title('ECG Lead 2');
    hold on; trigPoints = nan(length(waveform.ecg2),1);
    trigPoints(trigger.ecg2) = waveform.ecg2(trigger.ecg2);
    scatter(1:length(trigPoints),trigPoints,'filled','green');
    legend('ECG2','ECG2 Triggers'); hold off
    subplot(2,1,2); plot(waveform.ecg3); title('ECG Lead 3');
    hold on; trigPoints = nan(length(waveform.ecg3),1);
    trigPoints(trigger.ecg3) = waveform.ecg3(trigger.ecg3);
    scatter(1:length(trigPoints),trigPoints,'filled','green');
    legend('ECG3','ECG3 Triggers'); hold off
    SNRecg2 = snr(waveform.ecg2);
    SNRecg3 = snr(waveform.ecg3);
    if SNRecg2>-8 % this is an arbitray threshold value, can be changed
        fprintf('   ECG2: \t \t SIGNAL DETECTED\n');
    else
        fprintf('   ECG2: \t \t LIKELY NOISE\n');
    end 
    if SNRecg3>-8 % this is an arbitray threshold value, can be changed
        fprintf('   ECG3: \t \t SIGNAL DETECTED\n');
    else
        fprintf('   ECG3: \t \t LIKELY NOISE\n');
    end
    savefig('Raw_ecg_waveforms_triggers.fig')
end 



%% Catch
if isempty(filesLoaded) % if we didn't find anything to load..
    disp('No gating tracks were found');
    return % kick us out
end 
ecgData.prospectiveGating = prospectiveGating; % save prospectiveGating flag to compositive dataset
ecgData.encodeLength = encodeLength; % save encoding type to compositive dataset
ecgData.filesLoaded = filesLoaded; % save all files loaded into compositive dataset
save('ecgData.mat','ecgData'); % this file contains all loaded information
diary off % turn off diary, stop writing to ecgInformation.txt



%% Clear unnecessary variables
clear gate ans fid dirInfo acquisition_order ecg ecg_sorted fileNames
clear ecgIdx ecgName frame_order idx infoIdx infoName inter inter_index
clear inter_pos keepFiles looking name needINTER needNPROJ needORDER 
clear needTR nproj ORDER param1 param2 param3 param4 pos PR_ORDER
clear rhuser10 rhuser8 sub_order subproj time tline tr TR N sub_pos
clear yLimits timeres rr rrFull peaks peaksFull peaksECG rrECG
clear newLims minRR meanRR maxRR iterator i hFig dt trigPoints timeresFull
clear numMissed numMissedFull numEarly numEarlyFull missedHBidx 
clear expandedTime earlyTrigs earlyTrigIdx 
clear missedHBs gateOn respOn filtDir fullDir
clear SNRppg SNRecg2 SNRecg3 SNRresp badProspGating



%% Ancillary functions
function [RRs,peaks] = getRR(ecg)
    df = diff(ecg); % differentiate ecg with respect to time
    dt = mode(df); % get time resolution of ecg by finding most common value in derivative
    peaks = find(df<-12*dt); % find negative peaks, 12*dt is an arbitray value
    RRs = ecg(peaks); % get the actual rr value at the peak
end 

function [missedHBidx,earlyTrigIdx] = findMissedHB(rr,peaks)
    % Find missed heartbeats
    missedHBidx = []; % create missed hearbeat array
    earlyTrigIdx = []; % create early trigger array
    if length(rr)>100 % if we have more than 100 rr intervals
        for L = 1:length(rr)
            if L<25 % starting on the left side
                RR_window = linspace(1,50,50)'; % create sliding window over first 50 points
            elseif L>length(rr)-26 % in the middle, move the sliding window
                RR_window = linspace(length(rr)-50,length(rr)-1,50)'; % create sliding window of width 50
            else % if we're at the end of the rr vector
                RR_window = linspace(L-24,L+25,50)'; % create sliding window over last 50 points
            end
            sliding_RR = median(rr(RR_window)); % grab median rr from sliding window
            scale = rr(L)/sliding_RR; % find how the current rr interval compares to the median
            if scale > 1.5 % if scale is greater than 1.66 times the median rr, it is likely a missed heartbeat. 1.66 is arbitrary
                missedHBidx = [missedHBidx,peaks(L)]; % find index in the peak array and place in running array
            end 
            if scale <0.75 % if scale is less than 0.66 times the median rr, it is likely an early trigger. 0.66 is arbitrary
                earlyTrigIdx = [earlyTrigIdx,peaks(L)]; % find index in the peak array and place in running array
            end 
        end
    elseif length(rr)>40 % if we have less than 100 rr intervals but more than 40
        for L = 1:length(rr)
            if L<10
                RR_window = linspace(1,20,20)';
            elseif L>length(rr)-11
                RR_window = linspace(length(rr)-20,length(rr)-1,20)';
            else
                RR_window = linspace(L-9,L+10,20)';
            end
            sliding_RR = median(rr(RR_window));
            scale = rr(L)/sliding_RR;
            if scale > 1.66
                missedHBidx = [missedHBidx,peaks(L)];
            end 
            if scale <0.66
                earlyTrigIdx = [earlyTrigIdx,peaks(L)];
            end 
        end
    else % if we have a very small sample size (unreliable)
        for L = 1:length(rr)
            scale = rr(L)/median(rr);
            if scale > 1.66
                missedHBidx = [missedHBidx,peaks(L)];
            end 
            if scale <0.66
                earlyTrigIdx = [earlyTrigIdx,peaks(L)];
            end 
        end 
    end 
end 

function B = bitreverse(N)
    N_bits = ceil(log2(N));
    B_tmp = zeros(1, N);
    for j = 0 : N-1
         tmp = j;
         for i = 0 : N_bits
             B_tmp(j+1) = B_tmp(j+1) + 2^(N_bits - i - 1) * mod(tmp,2);
             tmp = floor(tmp/2);
         end
    end

    B(1) = B_tmp(1);
    curr_high = 0;

    for j = 2 : N
         tmp = N+1;
         diff_elt = 1;
         for i = 2 : N
             if B_tmp(i) > curr_high
                 diff = B_tmp(i) - curr_high;
                 if diff < tmp
                    diff_elt = i;
                    tmp = diff;
                 end
             end
         end
             curr_high = B_tmp(diff_elt);
             B(diff_elt) = j-1;
    end
end 