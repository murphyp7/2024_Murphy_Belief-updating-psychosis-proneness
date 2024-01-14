% Interpolate pupil time-series

function loop_interpSS(subj)

loadpath = '/mnt/homes/home024/pmurphy/Surprise_scz/pupil/3.matConverted/';
savepath = '/mnt/homes/home024/pmurphy/Surprise_scz/pupil/4.interpolated/';
addpath '/mnt/homes/home024/pmurphy/Surprise_scz/pupil/'

% get list of all subject IDs
files = dir([loadpath,'*.mat']);
for f = 1:length(files), allsubj{f}=files(f).name(2:4); end
allsubj = unique(allsubj);

% pull files for only current subject
subj = allsubj{subj};
files = dir([loadpath,'t',subj,'*.mat']);

% Loop through individual files
errfiles={};
for f = 1:length(files)
    fprintf('File %s...\n',files(f).name)
    inFile = [loadpath,files(f).name];
    outFile = [savepath,files(f).name];
    
    load(inFile)
    
    try
        % Run interpolation routine & add results to data structure
        [newpupil, newXgaze, newYgaze, newblinksmp, badsmp, params, h, bestchan] = blink_interpolate_1chan(data,files(f).name);
        data.pupil = newpupil;
        data.Xgaze = newXgaze;
        data.Ygaze = newYgaze;
        data.newblinksmp = newblinksmp;
        data.badsmp = badsmp;
        data.interp_params = params;
        data.bestchan = bestchan;
        
        data = rmfield(data,{'pupilL','pupilR'});
        
        save(outFile,'data')
        savefig(h, [savepath,'plots',filesep,files(f).name(1:end-4),'.fig'],'compact')
        close all
    catch
        errfiles{end+1,1} = files(f).name;
    end
end

if ~isempty(errfiles)
    save([savepath,'bad_files',filesep,subj,'_bad_files.mat'],'errfiles')
end

end