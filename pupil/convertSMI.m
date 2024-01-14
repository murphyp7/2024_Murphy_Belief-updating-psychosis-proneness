

% Convert from SMI .idf-converted text files to .asc files
loadpath = '/mnt/homes/home024/pmurphy/Surprise_scz/pupil/1.converted/';
savepath = '/mnt/homes/home024/pmurphy/Surprise_scz/pupil/2.ascConverted/';

addpath '/mnt/homes/home024/pmurphy/Surprise_scz/pupil/'
addpath '/mnt/homes/home024/pmurphy/Toolboxes/fieldtrip-20160221'  % tell Matlab where FieldTrip is
ft_defaults

% get files
files = dir([loadpath,'*.txt']);

% Loop through individual files
errfiles={};
for f = 1:length(files)
    fprintf('File %s...\n',files(f).name)
    matFile = [savepath,files(f).name(1:end-12),'.mat'];
    % Read the asc file into matlab & save
    if ~exist(matFile)  % ignores if file exists already
        try
            tic
            asc = read_SMI_asc(files(f).name);
            toc
            save(matFile,'asc')
        catch
            errfiles{end+1,1} = files(f).name;
        end
    end
end



