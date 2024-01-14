% Convert messy edf2asc conversions into my own format

clear, close all

loadpath = '/mnt/homes/home024/pmurphy/Surprise_scz/pupil/2.ascConverted/';
savepath = '/mnt/homes/home024/pmurphy/Surprise_scz/pupil/3.matConverted/';

addpath '/mnt/homes/home024/pmurphy/Surprise_scz/pupil/'
addpath '/mnt/homes/home024/pmurphy/Toolboxes/fieldtrip-20160221'  % tell Matlab where FieldTrip is  % 'tokenize' FT function is used

ft_defaults

% get files
files = dir([loadpath,'*.mat']);

% Loop through individual files
errfiles={}; trial_counts=[];
for f = 1%:length(files)
    fprintf('File %s...\n',files(f).name)
    inFile = [loadpath,files(f).name];
    outFile = [savepath,files(f).name];
    
    if ~exist(outFile)  % ignores if file exists already
        load(inFile)
        try
            
            tic
            % Create data and events structures, and matrices of blink/saccade times
            data = asc2dat(asc);
            
            % Manually adjust problem files
            if strcmp(files(f).name,'t021_2_5.mat')
                data.event = data.event(574:end);  % weird case where there's a bunch of events from different task at start of recording
            end
            
            % Create refined event structure
            data = surprise_trialfun(data,str2double(files(f).name(6)),str2double(files(f).name(8:end-4)));
            
            % Log [ID sess block ntrials]
            trial_counts(end+1,1:4) = [str2double(files(f).name(2:4)) ...
                str2double(files(f).name(6)) ...
                str2double(files(f).name(8:end-4)) ...
                size(data.event,1)];
            
            % Save output
            save(outFile,'data')
            
            toc
        catch
            errfiles{end+1,1} = files(f).name;
        end
    end
end

