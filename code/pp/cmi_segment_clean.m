function cmi_segment_clean(datafile, chan_info_file, run_on_server)
% cmi_segment_clean - grab raw data and cut out the segment with the data
% of interest for that task
%
%   datafile = full filename to the data of interest
%
%   Example:
%
%   datafile = '/Users/mlombardo/Dropbox/cmi_test/data/raw/NDARZC499NVX/RestingState/NDARZC499NVX_RestingState.mat';
%   chan_info_file = '/Users/mlombardo/Dropbox/cmi_test/code/GSN-HydroCel-129.sfp';
%   run_on_server = true;
%   cmi_segment_clean(datafile, chan_info_file, run_on_server)
%
%   datafile = '/media/DATA/RAW/cmihbn/data/raw/NDARZC499NVX/RestingState/NDARZC499NVX_RestingState.mat'
%   chan_info_file = '/media/DATA/RAW/cmihbn/code/GSN-HydroCel-129.sfp';
%   run_on_server = true;
%   cmi_segment_clean(datafile, chan_info_file, run_on_server)
%

%% ========================================================================
% load EEGlab
if run_on_server
    eeglabdir = '/home/mlombardo/eeglab_20201226';
elseif ~run_on_server
    eeglabdir = '/Users/mlombardo/Dropbox/matlab/eeglab_20201226';
end

rootpath = pwd;
cd(eeglabdir);
eeglab nogui;
cd(rootpath);



%% ========================================================================
% work out directories
disp(sprintf('Working on %s',datafile)); 
[fpath,fname,fext] = fileparts(datafile);
[subroot,task_name,fext] = fileparts(fpath);
[rawroot,subid,fext] = fileparts(subroot);
tmp_dir = fullfile(rawroot,'tmp','problematic_data');

%% ========================================================================
% run main functions to clean up data
try
    
    disp('Running clean up');
    if strcmp(task_name,'RestingState')
        EEG = clean_resting_state(datafile, chan_info_file);
    elseif strcmp(task_name,'Despicable')
        EEG = clean_despicable(datafile, chan_info_file);
    elseif strcmp(task_name,'Fractals')
        EEG = clean_fractals(datafile, chan_info_file);
    elseif strcmp(task_name,'Wimpy')
        EEG = clean_wimpy(datafile, chan_info_file);
    elseif strcmp(task_name,'Present')
        EEG = clean_present(datafile, chan_info_file);
    end % if strcmp(task_name,'RestingState')
    
    
    %% ====================================================================
    % save data to .set
    disp('Saving data');
    fname2save = fullfile(fpath,sprintf('%s.set',fname));
    pop_saveset(EEG, 'filename', fname2save);

catch
    
    % save an error message file to the directory to label it as
    % problematic
    error_message = '!!! Problem with data. Could not clean. !!!';
    disp(error_message);
    error_file2save = fullfile(fpath, 'cleaning_error.txt');
    unix_str = sprintf('echo %s >> %s', error_message, error_file2save);
    unix(unix_str);
    
    % move problematic data to the tmp directory
    tmp_sub_dir = fullfile(tmp_dir,subid);
    unix_str = sprintf('mkdir -p %s',tmp_sub_dir);
    unix(unix_str);
    unix_str = sprintf('mv %s %s',fpath,tmp_sub_dir);
    unix(unix_str);
    
    % check if subject directory is empty and if so, delete it
    dir_is_empty = check_dir(subroot);
    if dir_is_empty
        unix(sprintf('rm -Rf %s',subroot));
    end
    
end % try

end % function cmi_segment_clean


%% ========================================================================
% function check_dir to check if a directory is empty
function result = check_dir(dir2check)

d = dir(dir2check);
stuff2remove = {'.','..','.DS_Store'};
idx2remove = [];
for i = 1:length(d)
    if ismember(d(i).name,stuff2remove)
        idx2remove = [idx2remove,i];
    end
end % for i
d(idx2remove) = [];
result = isempty(d);

end % function add_channel_info



%% ========================================================================
% function to add channel info
function eeg_raw = add_channel_info(eeg_raw, chan_info_file)

disp('Adding channel info');
eeg_raw = pop_chanedit(eeg_raw, 'load',{chan_info_file,'filetype','sfp'});

end % function add_channel_info


%% ========================================================================
% function to add latencies to EEG.event field
function eeg_raw = add_latencies(eeg_raw)

disp('Adding latencies');

% sampling rate
sampling_rate = eeg_raw.srate;

% add latency and latency_sec fields to eeg_raw.event(i_event)
n_event = length(eeg_raw.event);
event_cell = cell(n_event,1);

for i_event = 1:n_event
    
    event_cell{i_event,1} = eeg_raw.event(i_event).type;
    
    % latency
    eeg_raw.event(i_event).latency = eeg_raw.event(i_event).sample;
    
    % latency_sec
    eeg_raw.event(i_event).latency_sec = (eeg_raw.event(i_event).sample -1 )/ sampling_rate;

end % for i_event = 1:n_event

% make urevent field
eeg_raw = eeg_checkset(eeg_raw, 'makeur');

end % function add_latencies


%% ========================================================================
function eeg_raw = clean_resting_state(datafile, chan_info_file)

% load data ---------------------------------------------------------------
load(datafile);
eeg_raw = EEG; % save an actual version of the raw loaded data before changing things inside the EEG structure

% -------------------------------------------------------------------------
% add channel info
eeg_raw = add_channel_info(eeg_raw, chan_info_file);

% -------------------------------------------------------------------------
% add latencies
eeg_raw = add_latencies(eeg_raw);

% -------------------------------------------------------------------------
% filter out junk

disp('Finding events of interest');

% event codes of interest (eoi) 
eoi = {'20  ', '30  '};

% grab event types
events2use = eeg_raw.event;
for ievent = 1:length(events2use)
    event_types{ievent} = events2use(ievent).type;
end % for ievent

% find events of interest (eoi)
mask = ismember(event_types, eoi);
eoi_idx = find(mask);

% find timepoint range of data to keep
starting_sample = eeg_raw.event(eoi_idx(1)).sample-1;
ending_sample = eeg_raw.event(eoi_idx(end)).sample;
timepoint_range = [starting_sample, ending_sample];
eeg_raw = pop_select(eeg_raw, 'point', timepoint_range);

end % function clean_resting_state


%% ========================================================================
% function eeg_raw = clean_despicable(datafile, chan_info_file)
function eeg_raw = clean_wimpy(datafile, chan_info_file)

ground_truth = 58699;
fudge_factor = 1000;

% load data ---------------------------------------------------------------
load(datafile);
eeg_raw = EEG; % save an actual version of the raw loaded data before changing things inside the EEG structure

% -------------------------------------------------------------------------
% add channel info
eeg_raw = add_channel_info(eeg_raw, chan_info_file);

% -------------------------------------------------------------------------
% add latencies
eeg_raw = add_latencies(eeg_raw);

% -------------------------------------------------------------------------
% filter out junk

disp('Finding events of interest');

% event codes of interest (eoi) 
eoi = {'81  ','82  ','83  ','84  ','101 ','102 ','103 ','104 '};
start_eoi = {'81  ','82  ','83  ','84  '};
end_eoi = {'101 ','102 ','103 ','104 '};

% grab event types
events2use = eeg_raw.event;
for ievent = 1:length(events2use)
    event_types{ievent} = events2use(ievent).type;
end % for ievent

% find events of interest (eoi)
mask = ismember(event_types, eoi);
eoi_idx = find(mask);
start_mask = ismember(event_types, start_eoi);
start_idx = find(start_mask);
end_mask = ismember(event_types, end_eoi);
end_idx = find(end_mask);
for i = 1:length(start_idx)
    duration = eeg_raw.event(start_idx(i)+1).sample - eeg_raw.event(start_idx(i)).sample;
    criteria2pass = ismember(duration,[(ground_truth-fudge_factor):(ground_truth+fudge_factor)]) & end_mask(start_idx(i)+1);
    if criteria2pass
        mask2use = zeros(1,length(event_types));
        mask2use(start_idx(i)) = 1;
        mask2use(start_idx(i)+1) = 1;
        mask2use = logical(mask2use);
        break
    end % if criteria2pass
end % for i

eoi_idx = find(mask2use);

% find timepoint range of data to keep
starting_sample = eeg_raw.event(eoi_idx(1)).sample-1;
ending_sample = eeg_raw.event(eoi_idx(end)).sample;
timepoint_range = [starting_sample, ending_sample];
eeg_raw = pop_select(eeg_raw, 'point', timepoint_range);

% end % function clean_despicable
end % function clean_wimpy


%% ========================================================================
% function eeg_raw = clean_fractals(datafile, chan_info_file)
function eeg_raw = clean_despicable(datafile, chan_info_file)

ground_truth = 85273;
fudge_factor = 1000;

% load data ---------------------------------------------------------------
load(datafile);
eeg_raw = EEG; % save an actual version of the raw loaded data before changing things inside the EEG structure

% -------------------------------------------------------------------------
% add channel info
eeg_raw = add_channel_info(eeg_raw, chan_info_file);

% -------------------------------------------------------------------------
% add latencies
eeg_raw = add_latencies(eeg_raw);

% -------------------------------------------------------------------------
% filter out junk

disp('Finding events of interest');

% event codes of interest (eoi) 
eoi = {'81  ','82  ','83  ','84  ','101 ','102 ','103 ','104 '};
start_eoi = {'81  ','82  ','83  ','84  '};
end_eoi = {'101 ','102 ','103 ','104 '};

% grab event types
events2use = eeg_raw.event;
for ievent = 1:length(events2use)
    event_types{ievent} = events2use(ievent).type;
end % for ievent

% find events of interest (eoi)
mask = ismember(event_types, eoi);
eoi_idx = find(mask);
start_mask = ismember(event_types, start_eoi);
start_idx = find(start_mask);
end_mask = ismember(event_types, end_eoi);
end_idx = find(end_mask);
for i = 1:length(start_idx)
    duration = eeg_raw.event(start_idx(i)+1).sample - eeg_raw.event(start_idx(i)).sample;
    criteria2pass = ismember(duration,[(ground_truth-fudge_factor):(ground_truth+fudge_factor)]) & end_mask(start_idx(i)+1);
    if criteria2pass
        mask2use = zeros(1,length(event_types));
        mask2use(start_idx(i)) = 1;
        mask2use(start_idx(i)+1) = 1;
        mask2use = logical(mask2use);
        break
    end % if criteria2pass
end % for i

eoi_idx = find(mask2use);

% find timepoint range of data to keep
starting_sample = eeg_raw.event(eoi_idx(1)).sample-1;
ending_sample = eeg_raw.event(eoi_idx(end)).sample;
timepoint_range = [starting_sample, ending_sample];
eeg_raw = pop_select(eeg_raw, 'point', timepoint_range);

% end % function clean_fractals
end % function clean_despicable


%% ========================================================================
% function eeg_raw = clean_wimpy(datafile, chan_info_file)
function eeg_raw = clean_fractals(datafile, chan_info_file)

ground_truth = 81493;
fudge_factor = 1000;

% load data ---------------------------------------------------------------
load(datafile);
eeg_raw = EEG; % save an actual version of the raw loaded data before changing things inside the EEG structure

% -------------------------------------------------------------------------
% add channel info
eeg_raw = add_channel_info(eeg_raw, chan_info_file);

% -------------------------------------------------------------------------
% add latencies
eeg_raw = add_latencies(eeg_raw);

% -------------------------------------------------------------------------
% filter out junk

disp('Finding events of interest');

% event codes of interest (eoi) 
eoi = {'81  ','82  ','83  ','84  ','101 ','102 ','103 ','104 '};
start_eoi = {'81  ','82  ','83  ','84  '};
end_eoi = {'101 ','102 ','103 ','104 '};

% grab event types
events2use = eeg_raw.event;
for ievent = 1:length(events2use)
    event_types{ievent} = events2use(ievent).type;
end % for ievent

% find events of interest (eoi)
mask = ismember(event_types, eoi);
eoi_idx = find(mask);
start_mask = ismember(event_types, start_eoi);
start_idx = find(start_mask);
end_mask = ismember(event_types, end_eoi);
end_idx = find(end_mask);
for i = 1:length(start_idx)
    duration = eeg_raw.event(start_idx(i)+1).sample - eeg_raw.event(start_idx(i)).sample;
    criteria2pass = ismember(duration,[(ground_truth-fudge_factor):(ground_truth+fudge_factor)]) & end_mask(start_idx(i)+1);
    if criteria2pass
        mask2use = zeros(1,length(event_types));
        mask2use(start_idx(i)) = 1;
        mask2use(start_idx(i)+1) = 1;
        mask2use = logical(mask2use);
        break
    end % if criteria2pass
end % for i

eoi_idx = find(mask2use);

% find timepoint range of data to keep
starting_sample = eeg_raw.event(eoi_idx(1)).sample-1;
ending_sample = eeg_raw.event(eoi_idx(end)).sample;
timepoint_range = [starting_sample, ending_sample];
eeg_raw = pop_select(eeg_raw, 'point', timepoint_range);

% end % function clean_wimpy
end % function clean_fractals



%% ========================================================================
function eeg_raw = clean_present(datafile, chan_info_file)

ground_truth = 101537;
fudge_factor = 1000;

% load data ---------------------------------------------------------------
load(datafile);
eeg_raw = EEG; % save an actual version of the raw loaded data before changing things inside the EEG structure

% -------------------------------------------------------------------------
% add channel info
eeg_raw = add_channel_info(eeg_raw, chan_info_file);

% -------------------------------------------------------------------------
% add latencies
eeg_raw = add_latencies(eeg_raw);

% -------------------------------------------------------------------------
% filter out junk

disp('Finding events of interest');

% event codes of interest (eoi) 
eoi = {'81  ','82  ','83  ','84  ','101 ','102 ','103 ','104 '};
start_eoi = {'81  ','82  ','83  ','84  '};
end_eoi = {'101 ','102 ','103 ','104 '};

% grab event types
events2use = eeg_raw.event;
for ievent = 1:length(events2use)
    event_types{ievent} = events2use(ievent).type;
end % for ievent

% find events of interest (eoi)
mask = ismember(event_types, eoi);
eoi_idx = find(mask);
start_mask = ismember(event_types, start_eoi);
start_idx = find(start_mask);
end_mask = ismember(event_types, end_eoi);
end_idx = find(end_mask);
for i = 1:length(start_idx)
    duration = eeg_raw.event(start_idx(i)+1).sample - eeg_raw.event(start_idx(i)).sample;
    criteria2pass = ismember(duration,[(ground_truth-fudge_factor):(ground_truth+fudge_factor)]) & end_mask(start_idx(i)+1);
    if criteria2pass
        mask2use = zeros(1,length(event_types));
        mask2use(start_idx(i)) = 1;
        mask2use(start_idx(i)+1) = 1;
        mask2use = logical(mask2use);
        break
    end % if criteria2pass
end % for i

eoi_idx = find(mask2use);

% find timepoint range of data to keep
starting_sample = eeg_raw.event(eoi_idx(1)).sample-1;
ending_sample = eeg_raw.event(eoi_idx(end)).sample;
timepoint_range = [starting_sample, ending_sample];
eeg_raw = pop_select(eeg_raw, 'point', timepoint_range);

end % function clean_present
