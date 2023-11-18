function eeg_pp = cmi_run_preproc(datafile, run_on_server)
% cmi_run_preproc - preprocess CMI EEG data for RestingState and Videos
%
%   INPUT
%       datafile = full filename to raw data *.set file
%       run_on_server = true or false, which tells the script whether
%                       you're running this on the lab's server or not
%
%   Example usage:
%
%       % Run on laptop
%       datafile = '/Users/mlombardo/Dropbox/cmi_test/data/raw/NDARET653TAM/RestingState/NDARET653TAM_RestingState.set';
%       run_on_server = false;
%       eeg_pp = cmi_run_preproc(datafile);
%
%       % Run on server
%       datafile = '/media/DATA/RAW/cmihbn/data/raw/NDARET653TAM/RestingState/NDARET653TAM_RestingState.set';
%       run_on_server = true;
%       eeg_pp = cmi_run_preproc(datafile);
%


%% initialize cfg structure
if run_on_server
    eeglab_dir = '/home/mlombardo/eeglab_20201226';
else
    eeglab_dir = '/Users/mlombardo/Dropbox/matlab/eeglab_20201226';
end % if run_on_server

cfg = initialize_cfg(datafile, run_on_server, eeglab_dir);

%% Run preprocessing pipeline with cmi_preproc_land
eeg_pp = cmi_preproc_land(cfg);

end % function cmi_run_preproc

%% Function to initialize cfg
function cfg = initialize_cfg(datafile, run_on_server, eeglab_dir)

% figure out paths and fname
[task_root, fname, fext] = fileparts(datafile);
[sub_root, dataType, fext] = fileparts(task_root);
[raw_root, subid, fext] = fileparts(sub_root);

% setup cfg structure
cfg.subid = subid;
cfg.dataType = dataType;
cfg.setfilename = datafile;
cfg.do_server = run_on_server;

% specifies where project directories, depending on whether they are on the server or not 
if cfg.do_server
    cfg.project_dir = '/media/DATA/RAW/cmihbn';
elseif ~cfg.do_server
    cfg.project_dir = '/Users/mlombardo/Dropbox/cmi_test';
end % if cfg.do_server

% other directories
cfg.raw_data_dir = fullfile(cfg.project_dir,'data','raw');
cfg.preproc_data_dir = fullfile(cfg.project_dir,'data','preproc');
cfg.code_dir = fullfile(cfg.project_dir,'code');
cfg.do_save_eeg_set = true;

% install EEGlab
cfg.eeglab_dir = eeglab_dir;
cd(cfg.eeglab_dir);
eeglab('nogui');
cd(cfg.code_dir);

% set up preprocessing parameters
pp.downsample_rate = 250; % sampling rate to downsample to in Hz

pp.hpf_cutoff = 1; % High pass filter in Hz
pp.lpf_cutoff = 80; % low pass filter cutoff in Hz

pp.line_noise_freq = 60; % line noise frequency in Hz
pp.line_noise_boundaries = 2; % how many Hz plus or minus the actual notch frequency to also include in the notch filter

% eye channels
pp.eye_channels = {'E128', 'E32', 'E25', 'E21', 'E127', 'E17', ...
    'E126', 'E14', 'E8', 'E1', 'E125'};

% reject channels
pp.do_chan_rejection = true;

% channels to reject
pp.chan_toreject = {
                'E127','E126',...
    'E25','E21','E17','E14','E8',...
    'E128','E32',          'E1','E125',...
    'E48','E43',            'E120','E119',...
    'E49',                          'E113',...
    'E56','E63',             'E99', 'E107',...
    'E68', 'E73', 'E81', 'E88', 'E94',...
               };

% interpolate channels
pp.do_chan_interp_pruning = true;

% channels to use for interpolation and then to prune before ICA
pp.chan_interp_prune =  {
    'E18','E10',...
    'E38','E121', ...
    'E44','E114', ...
    'E28','E117',...
    'E47','E98',...
    };

pp.iclabels_params = [0 0.2;0.7 1;0.7 1;NaN NaN;0.7 1;0.7 1;0.7 1]; % parameters for IClabels classification

% ASR parameters
pp.asr.iterations = 5;                                  % number of iterations to run ASR
pp.asr.chancrit = 0.8;                                  % 'ChannelCriterion'
pp.asr.lncrit = 4;                                      % 'LineNoiseCriterion'
pp.asr.flcrit = 5;                                      % 'FlatlineCriterion'
pp.asr.burstcrit = 20;                                  % 'BurstCriterion'
pp.asr.windowcrit = 'off';                              % 'WindowCriterion'
pp.asr.highpass = 'off';                                % 'Highpass'
pp.asr.burstrejection = 'on';                           % 'BurstRejection'
pp.asr.niter_bad_chan = round(pp.asr.iterations * 0.7); % n iterations for bad channels threshold

pp.do_plot_chan = false;
pp.do_plot_PSD = true;

pp.do_scorepoch = false;
pp.do_save_score = true;

% file stem for final preprocessed and denoised dataset 
pp.final_pp_dn_fstem = 'preproc_icaDenoised';
% file stem for final preprocessed but not denoised dataset
pp.final_pp_nodn_fstem = 'preproc_notDenoised';

% probability threshold for ICA component being classified as brain
pp.brain_prob_thresh = pp.iclabels_params(1,2);

% add pp to cfg
cfg.pp = pp;

% make subject's preproc directory
sub_pp_dir = fullfile(cfg.preproc_data_dir, cfg.subid, cfg.dataType);
unix_str = sprintf('mkdir -p %s', sub_pp_dir);
unix(unix_str);

end % function initialize_cfg


%%
function parameter_dump4report(cfg)

fprintf('Subject ID: %s \n',cfg.subid);
fprintf('Task: %s \n',cfg.dataType);
fprintf('File: %s \n',cfg.setfilename);
fprintf('\n');
fprintf('Downsampling rate: %d \n',cfg.pp.downsample_rate);
fprintf('High pass cutoff: %d \n',cfg.pp.hpf_cutoff);
fprintf('Low pass cutoff: %d \n',cfg.pp.lpf_cutoff);
fprintf('Line noise notch frequency: %d Hz + or - %d Hz \n', ...
    cfg.pp.line_noise_freq, ...
    cfg.pp.line_noise_boundaries);
fprintf('Outer channels to remove: \n');
cfg.pp.chan_toreject

if cfg.pp.do_chan_interp_pruning
    fprintf('Interpolate bad channels: true \n')
else
    fprintf('Interpolate bad channels: false \n')
end % if cfg.pp.do_chan_interp_pruning

fprintf('Channels to use for interpolation and then prune before ICA: \n');
cfg.pp.chan_interp_prune

fprintf('IClabels parameters to use when classifying components \n');
cfg.pp.iclabels_params

fprintf('ASR parameters: \n')
cfg.pp.asr

end % parameter_dump4report

%%
function [EEG] = cmi_preproc_land(cfg)
%
% CMI_PREPROC_LAND preprocess CMI data with LAND preprocessing pipeline
%

%% LOAD dataset
fprintf('Loading dataset \n')
eeg_struct = pop_loadset('filename',cfg.setfilename);
eeg_raw = eeg_struct;

subid = cfg.subid;
dataType = cfg.dataType;

%% Checks
sample_rate = eeg_struct.srate;
n_sample = eeg_struct.pnts;
n_chan = eeg_struct.nbchan;
n_chan_max = sqrt(n_sample/20);
% if n_sample > n_chan^2 * 20
%     disp([ num2str ' channels can be given as input to ICA'])
% else
%     sprintf('Number of channels for ICA should be reduced to %d' n_chan_max);
% end

%% Subset of 11 CHANNELS as EOG
chan_eye = cfg.pp.eye_channels;
eye_struct = pop_select(eeg_struct, 'channel', chan_eye);

%% Channel Removal
fprintf('Channel removal \n')
chan_toreject = cfg.pp.chan_toreject;
eeg_struct = pop_select(eeg_struct, 'nochannel', chan_toreject);
eeg_raw_chanred = eeg_struct;

%% Downsampling
fprintf('Downsampling \n')
eeg_struct = pop_resample(eeg_struct, cfg.pp.downsample_rate);
% fix latency so they are integers
for i = 1:length(eeg_struct.event)
    eeg_struct.event(i).latency = ceil(eeg_struct.event(i).latency);
end % for i
eeg_down = eeg_struct;

%% Band-pass filtering
fprintf('Band-pass filtering \n')
eeg_struct = pop_eegfiltnew(eeg_struct, cfg.pp.hpf_cutoff, [], [],0,[],0);
eeg_hpf = eeg_struct;
eeg_struct = pop_eegfiltnew(eeg_struct, [], cfg.pp.lpf_cutoff, [],0,[],0);
eeg_lpf = eeg_struct;

%% Remove line noise
fprintf('Remove line noise with notch filter \n')
eeg_notch = pop_eegfiltnew(eeg_struct, ...
    'locutoff',cfg.pp.line_noise_freq-cfg.pp.line_noise_boundaries, ...
    'hicutoff',cfg.pp.line_noise_freq+cfg.pp.line_noise_boundaries, ...
    'revfilt',1,'plotfreqz',1);

%% Iterative of bad channel and sample detection with cleanraw and ASR
n_iter = cfg.pp.asr.iterations;
fprintf('Iterative of bad channel and sample detection \n')

for i_iter = 1:n_iter
    
    sprintf('Iteration %d',i_iter);
    
    [eeg_cleanraw_iter, HP, BUR, bad_chan_cleanraw] = clean_artifacts(eeg_notch, ...
        'ChannelCriterion',cfg.pp.asr.chancrit, ...
        'LineNoiseCriterion',cfg.pp.asr.lncrit, ...
        'FlatlineCriterion',cfg.pp.asr.flcrit, ...
        'BurstCriterion',cfg.pp.asr.burstcrit, ...
        'WindowCriterion',cfg.pp.asr.windowcrit, ...
        'Highpass',cfg.pp.asr.highpass, ...
        'BurstRejection',cfg.pp.asr.burstrejection); 
    
    bad_chan_table(i_iter,:) = ~eeg_cleanraw_iter.etc.clean_channel_mask;
    bad_sample_table(i_iter,:) = ~eeg_cleanraw_iter.etc.clean_sample_mask;
end % for i_iter = 1:n_iter

bad_chan_idx = find(sum(bad_chan_table,1) > cfg.pp.asr.niter_bad_chan);  % bad channel identified in more than n iterations
bad_sample_idx = find(sum(bad_sample_table,1) > cfg.pp.asr.niter_bad_chan);  % bad channel identified in more than n iterations
bad_sample_perc = sum(bad_sample_table,2) ./ eeg_notch.pnts *100;

% final attempt to repair the bad samples
[eeg_cleanraw, HP, BUR, bad_chan_cleanraw] = clean_artifacts(eeg_notch, ...
    'ChannelCriterion',cfg.pp.asr.chancrit, ...
    'LineNoiseCriterion',cfg.pp.asr.lncrit, ...
    'FlatlineCriterion',cfg.pp.asr.flcrit, ...
    'BurstCriterion',cfg.pp.asr.burstcrit, ...
    'WindowCriterion',cfg.pp.asr.windowcrit, ...
    'Highpass',cfg.pp.asr.highpass, ...
    'BurstRejection','off');  

bad_chan_label = {};
counter = 1;
for i_chan = 1:length(eeg_cleanraw.etc.clean_channel_mask)
    if ~eeg_cleanraw.etc.clean_channel_mask(i_chan)
        bad_chan_label{1,counter} = eeg_notch.chanlocs(i_chan).labels;
        counter = counter+1;
    end % if ~eeg_cleanraw.etc.clean_channel_mask(i_chan)
end % for i_chan = 1:length(eeg_cleanraw.etc.clean_channel_mask)

% check which portion of the data is still considered as artifactual (even after ASR)
[eeg_cleanraw_ASR2, HP, BUR, bad_chan_cleanraw] = clean_artifacts(eeg_cleanraw, ...
    'ChannelCriterion',cfg.pp.asr.chancrit, ...
    'LineNoiseCriterion',cfg.pp.asr.lncrit, ...
    'FlatlineCriterion',cfg.pp.asr.flcrit, ...
    'BurstCriterion',cfg.pp.asr.burstcrit, ...
    'WindowCriterion',cfg.pp.asr.windowcrit, ...
    'Highpass',cfg.pp.asr.highpass, ...
    'BurstRejection',cfg.pp.asr.burstrejection);  


%% export csv file with num_bad_chan and pct_bad_samples
data2write = [{subid, dataType}, ... 
    num2cell([length(bad_chan_idx), nanmean(bad_sample_perc)])];
tab2write = cell2table(data2write, ...
    'VariableNames',{'subid','task','num_bad_channels','pct_bad_samples'});
fname2save = fullfile(cfg.preproc_data_dir, subid, dataType, ...
    sprintf('%s_%s_preproc_badSamplesChannels.csv',subid,dataType));
writetable(tab2write, fname2save,'Delimiter',',');


%% Figure plotting bad channels and samples
figure; set(gcf,'position',[10,10,1800,1000])
subplot(231)
imagesc(bad_chan_table')
xlabel('iteration');
ylabel('CHANNEL')
title(sprintf('n Bad Channels = %d',length(bad_chan_idx)));

subplot(234);
topoplot(sum(bad_chan_table,1), ...
    eeg_notch.chanlocs, ...
    'electrodes','labelpoint','chaninfo', ...
    eeg_notch.chaninfo) ;
title(bad_chan_label)

subplot(2,3,2:3)
imagesc(bad_sample_table)
ylabel('iteration');
xlabel('(BAD) SAMPLE  (10k samples = 40 sec)')
title(sprintf('Percentage of Bad Samples = %f', mean(bad_sample_perc)));

bad_sample_ASR2_table = double(~eeg_cleanraw_ASR2.etc.clean_sample_mask);
subplot(3,3,8:9)
imagesc(bad_sample_ASR2_table);
colormap(gca,'gray')
title('BAD segment  after initial ASR repair')

fname2save = fullfile(cfg.preproc_data_dir, subid, dataType, ...
    sprintf('%s_%s_preproc_badSamplesChannels.jpg',subid,dataType));
print(gcf,fname2save, '-djpeg','-noui');

%% Channel interpolation
fprintf('Channel interpolation \n')
eeg_cleanraw_badchan_interp = pop_interp(eeg_cleanraw, ...
    eeg_struct.chanlocs, ...
    'spherical');

%% Average re-referencing
fprintf('Average re-referencing \n'); 
eeg_cleanraw_avgref = pop_reref(eeg_cleanraw_badchan_interp, []);

%% Channel pruning before ICA
fprintf('Channel pruning before ICA \n')
chan_interp_prune = cfg.pp.chan_interp_prune;
eeg_cleanraw_avgref = pop_select(eeg_cleanraw_avgref, ...
    'nochannel', chan_interp_prune);

%% ICA
fprintf('Running ICA \n')
eeg_cleanraw_avgref_ICA = pop_runica(eeg_cleanraw_avgref, ...
    'icatype', ...
    'runica', ...
    'extended',1, ...
    'interrupt','on');

%% IClabels for classifying ICA components 
fprintf('Identifiy brain and non-brain ICA components with IClabels \n')
eeg_cleanraw_avgref_ICA = pop_iclabel(eeg_cleanraw_avgref_ICA, 'default');
% eeg_cleanraw_avgref_ICA = pop_icflag(eeg_cleanraw_avgref_ICA, ...
%     [0 0.2;0.7 1;0.7 1;NaN NaN;0.7 1;0.7 1;0.7 1]);
eeg_cleanraw_avgref_ICA = pop_icflag(eeg_cleanraw_avgref_ICA, ...
    cfg.pp.iclabels_params);

% find the brain and non-brain classified ICA components
brain_ica_components = find(eeg_cleanraw_avgref_ICA.reject.gcompreject == 0);
non_brain_ica_components = find(eeg_cleanraw_avgref_ICA.reject.gcompreject == 1);
n_icacomponents2reject = length(non_brain_ica_components);
n_icacomponents2keep = length(brain_ica_components);

%% Plot all ICA components and their probability of being classified as brain
% figure for plotting ICA components and their probability of being classified as brain
close all;
brain_prob_thresh = cfg.pp.brain_prob_thresh;
n_components = size(eeg_cleanraw_avgref_ICA.etc.ic_classification.ICLabel.classifications,1);
figure; plot(eeg_cleanraw_avgref_ICA.etc.ic_classification.ICLabel.classifications(:,1) ,1:n_components);
hold on; plot([brain_prob_thresh,brain_prob_thresh],ylim);
xlim([0,1]);
grid on;
ylabel('ICA Component');
xlabel('Brain Probability');
title('ICA Component Brain Probability');
yLims = ylim;
xLims = xlim;
text(xLims(2)-0.375, yLims(2)-10,sprintf('Brain components = %d',n_icacomponents2keep));
text(xLims(2)-0.375, yLims(2)-5,sprintf('Non-brain components = %d',n_icacomponents2reject));

fname2save = fullfile(cfg.preproc_data_dir, subid, dataType, ...
    sprintf('%s_%s_preproc_icacomponent_brain_probabilities.jpg',subid,dataType));
print(gcf,fname2save,'-djpeg','-noui');

%% Figure showing all ICA components and their classifications 
close all;

pop_viewprops(eeg_cleanraw_avgref_ICA, 0, [1:35], [2 80], []);
title('Components 1-35');
fname2save = fullfile(cfg.preproc_data_dir, subid, dataType, ...
    sprintf('%s_%s_preproc_ica_1_35.jpg',subid,dataType));
print(gcf,fname2save, '-djpeg','-noui');

pop_viewprops(eeg_cleanraw_avgref_ICA, 0, [36:70], [2 80], []);
title('Components 36-70');
fname2save = fullfile(cfg.preproc_data_dir, subid, dataType, ...
    sprintf('%s_%s_preproc_ica_36_70.jpg',subid,dataType));
print(gcf,fname2save, '-djpeg','-noui');

pop_viewprops(eeg_cleanraw_avgref_ICA, 0, [71:n_components], [2 80], []);
title('Components 71-93');
fname2save = fullfile(cfg.preproc_data_dir, subid, dataType, ...
    sprintf('%s_%s_preproc_ica_71_%d.jpg',subid,dataType,n_components));
print(gcf,fname2save, '-djpeg','-noui');

%% Project out all of the non-brain ICA components
fprintf('Projecting out all non-brain ICA components with IClabels \n');
if n_icacomponents2keep~=0
    eeg_cleanraw_avgref_nobadICA = pop_subcomp(eeg_cleanraw_avgref_ICA, non_brain_ica_components);
    EEG = eeg_cleanraw_avgref_nobadICA;
else
    fprintf('No ICA components classified as brain! \n');
    fprintf('Final data will be data preprocessed with all steps up until ICA. \n');
    EEG = eeg_cleanraw_avgref_ICA;
end % if n_icacomponents2keep~=0

%% Plot figure showing PSDs before and after ICA denoising
close all;
figure; set(gcf,'color','white');
if n_icacomponents2keep~=0

    subplot(1,2,1);
    pop_spectopo(eeg_cleanraw_avgref_ICA, 1, [ ], ...
        'EEG' , 'percent', 50, 'freq', [8 13 20], ...
        'freqrange',[2 80],'electrodes','on');
    
    subplot(1,2,2);
    pop_spectopo(eeg_cleanraw_avgref_nobadICA, 1, [ ], ...
        'EEG' , 'percent', 50, 'freq', [8 13 20], ...
        'freqrange',[2 80],'electrodes','on');
    sgtitle('PSD Before and After ICA');
    
else
    pop_spectopo(eeg_cleanraw_avgref_ICA, 1, [ ], ...
        'EEG' , 'percent', 50, 'freq', [8 13 20], ...
        'freqrange',[2 80],'electrodes','on');
    title('PSD Before ICA');
end 

fname2save = fullfile(cfg.preproc_data_dir, subid, dataType, ...
    sprintf('%s_%s_preproc_psd.jpg',subid,dataType));
print(gcf,fname2save, '-djpeg','-noui');

%% Save final preproc data to disk
if n_icacomponents2keep~=0
    fname2save = fullfile(cfg.preproc_data_dir, subid, dataType, ...
        sprintf('%s_%s_%s.set',subid, dataType, cfg.pp.final_pp_dn_fstem));
    pop_saveset(eeg_cleanraw_avgref_nobadICA, 'filename', fname2save);
else
    fname2save = fullfile(cfg.preproc_data_dir, subid, dataType, ...
        sprintf('%s_%s_%s.set',subid, dataType, cfg.pp.final_pp_nodn_fstem));
    pop_saveset(eeg_cleanraw_avgref_ICA, 'filename', fname2save);
end

%% Parameter dump for report
fprintf('Dumping out all parameters used for preprocessing \n');
parameter_dump4report(cfg);

fprintf('... Done! \n');

end % function cmi_preproc_land
