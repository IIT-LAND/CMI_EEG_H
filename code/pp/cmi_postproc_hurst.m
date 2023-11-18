function [results, outfile] = cmi_postproc_hurst(datafile, output_to_csv)
% cmi_postproc_hurst - compute the Hurst exponent on a preprocessed *.set file
%
% INPUT
%   datafile = full filename to a preprocessed *.set file.
%   output_to_csv = set to true to write the output to a csv file. Otherwise set to false.
%
% OUTPUT
%   results = a [n_electrodes, n_blocks+1] table of H values. First column
%   are electrode names.
%
%   This function also writes out a *.csv file if output_to_csv is true.
%
% Example usage:
%
%   datafile = '/media/DATA/RAW/cmihbn/data/preproc/NDARAM848GTE/RestingState/NDARAM848GTE_RestingState_cleanraw_avgref_nobadICA.set';
%   output_to_csv = true;
%   results = cmi_postproc_hurst(datafile, output_to_csv);
%
% -- written by nbertelsen and mvlombardo
%

%% Step 1: Get meta data, install EEGlab, and get H parameters
% get meta data
meta_data = cmi_get_meta_data(datafile);

% install EEGlab
cmi_install_eeglab(meta_data);

% get H parameters
Hparams = cmi_hurst_parameters;


%% Step 2: Prepare data
message2print = sprintf('Computing H on %s %s ',meta_data.sub_name, meta_data.task_name);
disp(message2print)

% prepare data
[data_segments, block_info] = cmi_prepare_data(meta_data);


%% Step 3: compute H
block_names = block_info.block_names;
electrodes = data_segments(1).electrodes;
results = cmi_compute_H(data_segments, block_names, electrodes, Hparams);


%% Step 4: write data out to a file, if needed
if output_to_csv
    cmi_write_output(results, meta_data);
end % if output_to_csv

outfile = meta_data.out_file;

end % function cmi_postproc_hurst


%% function cmi_hurst_parameters
% declare parameters used in H computation
function results = cmi_hurst_parameters

results.filter = 'haar';
results.lb = [-0.5, 0];
results.ub = [1.5, 10];
results.verbose = false;

end % function cmi_hurst_parameters


%% function cmi_get_meta_data
% get meta data like path info, subjectId, task name, and output filename
function results = cmi_get_meta_data(datafile)

% parse out paths, subjectId, and task
% parse apart the input datafile
results.datafile = datafile;
[results.task_path, results.fname, results.set_ext] = fileparts(datafile);
% parse apart the task_path
[results.sub_path, results.task_name] = fileparts(results.task_path);
% parse apart the sub_path
[results.preproc_path, results.sub_name] = fileparts(results.sub_path);
% parse apart the sub_path
[results.data_path] = fileparts(results.preproc_path);
% parse apart the sub_path
[results.root_path] = fileparts(results.data_path);
results.code_path = fullfile(results.root_path,'code');
if strcmp(results.root_path,'/media/DATA/RAW/cmihbn')
  results.eeglab_path = '/home/mlombardo/eeglab_20201226';
else
  results.eeglab_path = fullfile(results.code_path,'toolboxes','eeglab_20201226');
end
results.nfmaster_dir = fullfile(results.code_path,'toolboxes','nonfractal-master');
results.wmtsa_dir = fullfile(results.code_path,'toolboxes','wmtsa-matlab-0.2.6');

% make out_path and out_file
results.out_path = fullfile(results.data_path, 'postproc', results.sub_name, results.task_name, 'H');
results.out_file = fullfile(results.out_path, sprintf('%s_%s_H.csv', results.sub_name, results.task_name));

end % function cmi_get_meta_data


%% function cmi_install_eeglab
% install eeglab
function cmi_install_eeglab(meta_data)

cd(meta_data.eeglab_path);
eeglab nogui;
cd(meta_data.code_path);

% add nonfractal and wmtsa toolboxes
addpath(genpath(meta_data.nfmaster_dir));
addpath(genpath(meta_data.wmtsa_dir));

end % function cmi_install_eeglab

%% function cmi_prepare_data
% load in data, check it, and cut into data segments
function [data_segments, block_info] = cmi_prepare_data(meta_data)

% load data
EEG = pop_loadset(meta_data.datafile);

% check to make sure data is double format because bfn_mfin_ml won't run properly otherwise
if ~strcmp(class(EEG.data),'double')
  EEG.data = double(EEG.data);
end % if ~strcmp(class(EEG.data),'double')

if strcmp(meta_data.task_name, 'RestingState')
    % cut data into eyes open and closed segments
    data_segments = cmi_cut_rs_segments(EEG);
else
    % grab data segment as entire time-series
    data_segments.label = meta_data.task_name;
    data_segments.data = EEG.data;
    data_segments.electrodes = {EEG.chanlocs.labels};
end % if strcmp(task_name, 'RestingState')

% get block_info
block_info.block_names = {data_segments.label};
block_info.n_blocks = length(block_info.block_names);

end % function cmi_prepare_data

%% function for cutting resting state into eyes open vs closed segements
% cut up resting state data into blocks of eyes open versus closed
function data_segments = cmi_cut_rs_segments(EEG)

% assumes you have 5 blocks of eyes open and eyes closed
n_blocks = 5;

% initialize structure
data_segments = struct([]);

% Eyes open
onsets = [2:2:10];
offsets = [3:2:11];
for iblock = 1:n_blocks
    data_segments(iblock).label = sprintf('open%d',iblock);
    samples2use = [ceil(EEG.event(onsets(iblock)).latency):(floor(EEG.event(offsets(iblock)).latency)-1)];
    data_segments(iblock).data = EEG.data(:,samples2use);
    data_segments(iblock).electrodes = {EEG.chanlocs.labels};
end % for i = 1:n_blocks

% Eyes closed
counter = length(data_segments);
onsets = [3:2:11];
offsets = [4:2:12];
for iblock = 1:n_blocks
    counter = counter+1;
    data_segments(counter).label = sprintf('closed%d',iblock);
    samples2use = [ceil(EEG.event(onsets(iblock)).latency):(floor(EEG.event(offsets(iblock)).latency)-1)];
    data_segments(counter).data = EEG.data(:,samples2use);
    data_segments(counter).electrodes = {EEG.chanlocs.labels};
end % for i = 1:n_blocks

end % function cmi_cut_rs_segments


%% function cmi_compute_H
% main code that computes H
function results = cmi_compute_H(data_segments, block_names, electrodes, Hparams)

n_blocks = length(block_names);
n_electrodes = length(electrodes);

% pre-allocate H
H = nan(n_electrodes, n_blocks);

for iblock = 1:n_blocks

    message2print = sprintf('Block %d ',iblock);
    disp(message2print)

    % grab data to use for H computation
    data2use = data_segments(iblock).data';

    %% compute Hurst exponent
    H(:, iblock) = bfn_mfin_ml(data2use, ...
        'filter', Hparams.filter, ...
        'lb', Hparams.lb, ...
        'ub', Hparams.ub, ...
        'verbose',Hparams.verbose);

end % for iblock = 1:n_blocks

% make output table
Htab = array2table(H, 'VariableNames',block_names);
results = [cell2table(electrodes','VariableNames',{'Electrodes'}), Htab];

end % function cmi_compute_H

%% function cmi_write_output
% write H results to a csv file
function cmi_write_output(results, meta_data)

% make the out_path and out_file
unix_str = sprintf('mkdir -p %s', meta_data.out_path);
unix(unix_str);

% write out_file
message2print = sprintf('Writing %s',meta_data.out_file);
disp(message2print)
writetable(results, meta_data.out_file);

end % function cmi_write_output
