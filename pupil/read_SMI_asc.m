function [asc] = read_SMI_asc(filename)

% READ_SMI_ASC reads the header information, input triggers, messages
% and all data points from an SMI converted *.asc file
%
% Adapted from Fieldtrip's READ_EYELINK_ASC with edits from NK, AU & PM.


fprintf('reading in %s ...\n', filename);
fid = fopen(filename, 'rt');

asc = struct;
asc.header  = {};
asc.msg     = {};
asc.datcols = [];
asc.fsample = [];
asc.dat     = [];

current=0;
tmpdat={};

while ~feof(fid)
    tline = fgetl(fid);  % read current line
            
    if regexp(tline, '##')  % store all HEADER INFO
        asc.header = cat(1, asc.header, {tline});
        
        if regexp(tline, 'Sample Rate:')  % pull sampling rate
            tok = tokenize(tline);
            asc.fsample = str2num(tok{4});
        end
        
    elseif regexp(tline, 'MSG')  % store all MESSAGES
        asc.msg = cat(1, asc.msg, {tline});
        
    elseif regexp(tline, 'Time')  % store all DATA COLUMN HEADERS
        tmpcols = tokenize(tline,char(9));  % char(9) is tab delimiter
        asc.datcols = tmpcols(~strcmp(tmpcols,'Type'));  % this excludes the 'Type' column and keep only data
        
    elseif regexp(tline, 'SMP')  % store all SAMPLES
        current = current+1;
        
        tmp  = tokenize(tline,char(9));  % read all sample rows into cell array (will include 'Type' column, to be removed later)
        
        if size(asc.dat,1)<current
            % increase the allocated number of samples - preallocating in
            % this way speeds up x2
            asc.dat(end+10000,1:length(tmp)-1) = 0;
        end
        
        asc.dat(current,:) = [cellfun(@str2double, tmp(~strcmp(tmpcols,'Type')))];  % convert string to double and add to data array        
%         tmpdat(end+1,1:length(tmp)) = tmp;
        
    else
        % all other lines are not parsed
    end
end

% close the file
fclose(fid);

% remove the samples that were not filled with real data
asc.dat = asc.dat(1:current,:);

% add data to asc.dat
% for col = 1:length(tmpcols)
%     if ~isnan(str2double(tmpdat{1,col}));  % this will exclude the 'Type' column and keep only data
%         asc.datcols{1,end+1} = tmpcols{col};
%         asc.dat(:,end+1) = cellfun(@str2double, tmpdat(:,col));
%     end
% end

end

