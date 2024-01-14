function data = asc2dat(asc)
% takes asc data from SMI file and converts this into events and data structure

% interpolate to account for missing samples/uneven sampling rate
t1 = asc.dat(1,strcmp('Time',asc.datcols));  % get first time point (at microsecond precision)
ts = t1:(1000000/asc.fsample):asc.dat(end,strcmp('Time',asc.datcols));  % define time-points for interpolation
tempdat=[];
for col = 1:length(asc.datcols)
    if strcmp('Time',asc.datcols{col})
        tempdat(:,col) = ts;
    else
        tempdat(:,col) = interp1(asc.dat(:,strcmp('Time',asc.datcols)),asc.dat(:,col),ts);
    end
end
asc.dat = tempdat;

% make data struct
data                = [];

% create event structure for messages
evcell = cell(length(asc.msg),1);
event = struct('type', evcell, 'sample', evcell, 'value', evcell, 'offset', evcell, 'duration', evcell );

for i=1:length(asc.msg)
    
    strtok = tokenize(asc.msg{i});
    event(i).type = strtok{6};
    
    % match the message to its sample
    smpstamp = dsearchn(asc.dat(:,1), str2double(strtok{1}));
    % find closest sample index of trigger in ascii dat
    
    if ~isempty(smpstamp)
        event(i).sample = smpstamp(1);
    else % if no exact sample was found
        warning('no sample found');
    end
    event(i).value = asc.msg{i};
end

data.event = event;  % add events into data structure

% add data into data structure
% important: match the right data chans to their corresponding labels...
data.Xgaze          = asc.dat(:,strcmp('L POR X [px]',asc.datcols))';  % taking estimated gaze position from left eye - should match exactly that for right eye
data.Ygaze          = asc.dat(:,strcmp('L POR Y [px]',asc.datcols))';

if sum(strcmp('L Mapped Diameter [mm]',asc.datcols))>0
    data.pupilL      = asc.dat(:,strcmp('L Mapped Diameter [mm]',asc.datcols))';
    data.pupilR      = asc.dat(:,strcmp('R Mapped Diameter [mm]',asc.datcols))';
else   % if mapped diameter doesn't exist (as is case for small number of files), take untransformed pixels instead
    data.pupilL      = asc.dat(:,strcmp('L Dia X [px]',asc.datcols))';
    data.pupilR      = asc.dat(:,strcmp('R Dia X [px]',asc.datcols))';
    data.noMapped    = 1;
end

data.fsample        = asc.fsample;
data.times          = 0:1/data.fsample:length(asc.dat(:,1))/data.fsample-1/data.fsample;
data.sampleinfo     = [1 length(asc.dat(:,1))];

if data.fsample ~= 250
    warning('pupil not sampled with 250Hz');
end

data.blinksmp = [];
data.saccsmp = [];

end