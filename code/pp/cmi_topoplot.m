function cmi_topoplot(datafile)
% cmi_topoplot - take a postprocessing metric (e.g., intsc, H) and make a
% topoplot and save it as a pdf
%
%   INPUT
%       datafile = full filename to the *.csv file you want to plot
%
%   Example usage:
%       datafile = '/media/DATA/RAW/cmihbn/data/postproc/NDARVB151GXF/Present/intsc/NDARVB151GXF_Present_intsc.csv
%       cmi_topoplot(datafile);
%
%   written by mvlombardo
%

% get meta data
meta_data = cmi_get_meta_data(datafile);

% install EEGlab
cmi_install_eeglab(meta_data);

% load chanlocs
try
  EEG = pop_loadset(meta_data.preproc_setfile1);
catch
  EEG = pop_loadset(meta_data.preproc_setfile2);
end
chanlocs = EEG.chanlocs;

% grab data to plot
data2use = readtable(meta_data.datafile);

% make the figure and save it
if strcmp(meta_data.metric_name,'fooof')
    cmi_make_topoplot_fooof(data2use, chanlocs, meta_data.task_name, meta_data.sub_name, meta_data.metric_name, meta_data.out_fig_file)
else
    cmi_make_topoplot_other(data2use, chanlocs, meta_data.task_name, meta_data.sub_name, meta_data.metric_name, meta_data.out_fig_file)
end % if strcmp
end % function cmi_topoplot


%% function cmi_get_meta_data
% get meta data like path info, subjectId, task name, and output filename
function results = cmi_get_meta_data(datafile)

% parse out paths, subjectId, and task
% parse apart the input datafile
results.datafile = datafile;
[results.out_path, results.fname, results.set_ext] = fileparts(datafile);
% get task_path
[results.task_path, results.metric_name] = fileparts(results.out_path);
[results.sub_path, results.task_name] = fileparts(results.task_path);

% parse apart the sub_path
[results.postproc_path, results.sub_name] = fileparts(results.sub_path);
% get data_path
[results.data_path] = fileparts(results.postproc_path);

% get root_path
[results.root_path] = fileparts(results.data_path);
results.code_path = fullfile(results.root_path,'code');
if strcmp(results.root_path,'/media/DATA/RAW/cmihbn')
  results.eeglab_path = '/home/mlombardo/eeglab_20201226';
else
  results.eeglab_path = fullfile(results.code_path,'toolboxes','eeglab_20201226');
end
% results.nfmaster_dir = fullfile(results.code_path,'toolboxes','nonfractal-master');
% results.wmtsa_dir = fullfile(results.code_path,'toolboxes','wmtsa-matlab-0.2.6');

results.preproc_path = fullfile(results.data_path,'preproc',results.sub_name, results.task_name);
% results.preproc_setfile = fullfile(results.preproc_path,sprintf('%s_%s_cleanraw_avgref_nobadICA.set',results.sub_name, results.task_name));
results.preproc_setfile1 = fullfile(results.preproc_path,sprintf('%s_%s_preproc_icaDenoised.set',results.sub_name, results.task_name));
results.preproc_setfile2 = fullfile(results.preproc_path,sprintf('%s_%s_preproc_notDenoised.set',results.sub_name, results.task_name));


% make out_fig_file
results.out_fig_file = fullfile(results.out_path, sprintf('%s_%s_%s_topoplot.pdf', results.sub_name, results.task_name, results.metric_name));

end % function cmi_get_meta_data


%% function cmi_install_eeglab
% install eeglab
function cmi_install_eeglab(meta_data)

cd(meta_data.eeglab_path);
eeglab nogui;
cd(meta_data.code_path);

end % function cmi_install_eeglab


%% function cmi_reduce_fooof
function data2use = cmi_def_freq_bands(data2use)

delta_power = [0,3.99];
theta_power = [4,7.99];
alpha_power = [8,11.99];
beta_power = [12,29.99];
gamma_power = [30,50];

data2use.band = repmat({'band'},size(data2use,1),1);
mask = data2use.peak_frequency>=delta_power(1) & data2use.peak_frequency<=delta_power(2);
data2use.band(mask) = {'delta'};
mask = data2use.peak_frequency>=theta_power(1) & data2use.peak_frequency<=theta_power(2);
data2use.band(mask) = {'theta'};
mask = data2use.peak_frequency>=alpha_power(1) & data2use.peak_frequency<=alpha_power(2);
data2use.band(mask) = {'alpha'};
mask = data2use.peak_frequency>=beta_power(1) & data2use.peak_frequency<=beta_power(2);
data2use.band(mask) = {'beta'};
mask = data2use.peak_frequency>=gamma_power(1) & data2use.peak_frequency<=gamma_power(2);
data2use.band(mask) = {'gamma'};

end % cmi_def_freq_bands

%% function cmi_reduce_fooof
function results = cmi_reduce_fooof(data2use)

electrodes = unique(data2use.electrode);
block_names = unique(data2use.block);

% get oscillatory frequency bands
data2use = cmi_def_freq_bands(data2use);

% pre-allocate
offset_res = nan(length(electrodes),length(block_names));
exponent_res = nan(length(electrodes),length(block_names));
delta_freq_res = nan(length(electrodes),length(block_names));
delta_amp_res = nan(length(electrodes),length(block_names));
delta_bandwidth_res = nan(length(electrodes),length(block_names));
theta_freq_res = nan(length(electrodes),length(block_names));
theta_amp_res = nan(length(electrodes),length(block_names));
theta_bandwidth_res = nan(length(electrodes),length(block_names));
beta_freq_res = nan(length(electrodes),length(block_names));
beta_amp_res = nan(length(electrodes),length(block_names));
beta_bandwidth_res = nan(length(electrodes),length(block_names));
alpha_freq_res = nan(length(electrodes),length(block_names));
alpha_amp_res = nan(length(electrodes),length(block_names));
alpha_bandwidth_res = nan(length(electrodes),length(block_names));
gamma_freq_res = nan(length(electrodes),length(block_names));
gamma_amp_res = nan(length(electrodes),length(block_names));
gamma_bandwidth_res = nan(length(electrodes),length(block_names));

for i = 1:length(electrodes)

    electrode = electrodes{i};
    electrode_mask = ismember(data2use.electrode,electrode);

    for iblock = 1:length(block_names)

        block = block_names{iblock};
        block_mask = ismember(data2use.block,block);
        mask = electrode_mask & block_mask;
        offset_res(i,iblock) = unique(data2use.aper_offset(mask));
        exponent_res(i,iblock) = unique(data2use.aper_exponent(mask));

        band_mask = ismember(data2use.band,'delta');
        mask2use = mask & band_mask;
        delta_freq_res(i,iblock) = nanmean(data2use.peak_frequency(mask2use));
        delta_amp_res(i,iblock) = nanmean(data2use.peak_amplitude(mask2use));
        delta_bandwidth_res(i,iblock) = nanmean(data2use.peak_bandwidth(mask2use));

        band_mask = ismember(data2use.band,'theta');
        mask2use = mask & band_mask;
        theta_freq_res(i,iblock) = nanmean(data2use.peak_frequency(mask2use));
        theta_amp_res(i,iblock) = nanmean(data2use.peak_amplitude(mask2use));
        theta_bandwidth_res(i,iblock) = nanmean(data2use.peak_bandwidth(mask2use));

        band_mask = ismember(data2use.band,'alpha');
        mask2use = mask & band_mask;
        alpha_freq_res(i,iblock) = nanmean(data2use.peak_frequency(mask2use));
        alpha_amp_res(i,iblock) = nanmean(data2use.peak_amplitude(mask2use));
        alpha_bandwidth_res(i,iblock) = nanmean(data2use.peak_bandwidth(mask2use));

        band_mask = ismember(data2use.band,'beta');
        mask2use = mask & band_mask;
        beta_freq_res(i,iblock) = nanmean(data2use.peak_frequency(mask2use));
        beta_amp_res(i,iblock) = nanmean(data2use.peak_amplitude(mask2use));
        beta_bandwidth_res(i,iblock) = nanmean(data2use.peak_bandwidth(mask2use));

        band_mask = ismember(data2use.band,'gamma');
        mask2use = mask & band_mask;
        gamma_freq_res(i,iblock) = nanmean(data2use.peak_frequency(mask2use));
        gamma_amp_res(i,iblock) = nanmean(data2use.peak_amplitude(mask2use));
        gamma_bandwidth_res(i,iblock) = nanmean(data2use.peak_bandwidth(mask2use));

    end % for iblock

end % for i

results.aperiodic_offset = array2table(offset_res,'VariableNames',block_names);
results.aperiodic_exponent = array2table(exponent_res,'VariableNames',block_names);
results.delta_freq = array2table(delta_freq_res,'VariableNames',block_names);
results.delta_amp = array2table(delta_amp_res,'VariableNames',block_names);
results.delta_bandwidth = array2table(delta_bandwidth_res,'VariableNames',block_names);
results.theta_freq = array2table(theta_freq_res,'VariableNames',block_names);
results.theta_amp = array2table(theta_amp_res,'VariableNames',block_names);
results.theta_bandwidth = array2table(theta_bandwidth_res,'VariableNames',block_names);
results.alpha_freq = array2table(alpha_freq_res,'VariableNames',block_names);
results.alpha_amp = array2table(alpha_amp_res,'VariableNames',block_names);
results.alpha_bandwidth = array2table(alpha_bandwidth_res,'VariableNames',block_names);
results.beta_freq = array2table(beta_freq_res,'VariableNames',block_names);
results.beta_amp = array2table(beta_amp_res,'VariableNames',block_names);
results.beta_bandwidth = array2table(beta_bandwidth_res,'VariableNames',block_names);
results.gamma_freq = array2table(gamma_freq_res,'VariableNames',block_names);
results.gamma_amp = array2table(gamma_amp_res,'VariableNames',block_names);
results.gamma_bandwidth = array2table(gamma_bandwidth_res,'VariableNames',block_names);

end % function cmi_reduce_fooof


%% function for cmi_fooof_plot
function cmi_fooof_plot(data2use, block_names, chanlocs, ...
    sub_name, task_name, metric_name, band_name, measure_name, fname2save)

figure; set(gcf,'color','white');
if strcmp(measure_name,'offset')
    colormap2use = lbmap(100,'RedBlue');
else
    colormap2use = flipud(lbmap(100,'RedBlue'));
end
colorlimits_percentiles = [10,90]; % percentiles for color scaling

if strcmp(task_name,'RestingState')

%     block_names = {'open1','open2','open3','open4','open5','closed1','closed2','closed3','closed4','closed5'};
    concat_data = [];
    for i = 1:length(block_names)
        concat_data = [concat_data; data2use.(block_names{i})];
    end % for i = 1:length(block_names)
    colorlimits = prctile(concat_data,colorlimits_percentiles);
    tlo = tiledlayout(2,length(block_names)/2);

    for i = 1:length(block_names)
        h(i) = nexttile(tlo);
        topoplot(data2use.(block_names{i}), chanlocs, 'colormap', colormap2use, 'maplimits', colorlimits,'plotrad',0.5);
        title(block_names{i});

        if i==length(block_names)
            cbh = colorbar(h(i));
        end % if i==10

    end % for i = 1:length(block_names)

    sgtitle(sprintf('%s %s %s %s %s',sub_name, task_name, metric_name, band_name, measure_name));

else

    data2plot = data2use.(task_name);
    colorlimits = prctile(data2plot,colorlimits_percentiles);
    topoplot(data2plot, chanlocs, 'colormap',colormap2use,'maplimits',colorlimits); colorbar;
    title(sprintf('%s %s %s  %s %s',sub_name, task_name, metric_name, band_name, measure_name));

end % if strcmp(task_name,'RestingState')

% save figure
saveas(gcf, fname2save);

end % function cmi_fooof_plot





%% function for making the topoplot
function cmi_make_topoplot_fooof(data2use, chanlocs, task_name, sub_name, metric_name, out_fig_file)

fooof_results = cmi_reduce_fooof(data2use);

if strcmp(task_name,'RestingState')
    block_names = {'open1','open2','open3','open4','open5','closed1','closed2','closed3','closed4','closed5'};
else
    block_names = task_name;
end % if strcmp(task_name,'RestingState')

% oscillations in different frequency bands ---------------------------
freq_bands = {'delta','theta','alpha','beta','gamma'};
measure_names = {'freq','amp','bandwidth'};

for iband = 1:length(freq_bands)
    band_name = freq_bands{iband};

    for imeas = 1:length(measure_names)
        measure_name = measure_names{imeas};
        tab2use = fooof_results.(sprintf('%s_%s',band_name, measure_name));

        out_path = fileparts(out_fig_file);
        fname2save = fullfile(out_path, sprintf('%s_%s_%s_%s_%s_topoplot.pdf', ...
            sub_name, task_name, metric_name, ...
            band_name,measure_name));

        figure; set(gcf,'color','white');
        colormap2use = flipud(lbmap(100,'RedBlue'));
        colorlimits_percentiles = [10,90]; % percentiles for color scaling

        if strcmp(task_name,'RestingState')

            block_names = {'open1','open2','open3','open4','open5','closed1','closed2','closed3','closed4','closed5'};
            concat_data = [];
            for i = 1:length(block_names)
                concat_data = [concat_data; tab2use.(block_names{i})];
            end % for i = 1:length(block_names)
            colorlimits = prctile(concat_data,colorlimits_percentiles);
            tlo = tiledlayout(2,length(block_names)/2);

            for i = 1:length(block_names)
                h(i) = nexttile(tlo);
                if (sum(isnan(tab2use.(block_names{i})))==length(tab2use.(block_names{i})))
                    data2plot = repmat(0, length(tab2use.(block_names{i})),1);
                else
                    data2plot = tab2use.(block_names{i});
                end
                topoplot(data2plot, chanlocs, 'colormap', colormap2use, 'maplimits', colorlimits);
                title(block_names{i});

                if i==length(block_names)
                    cbh = colorbar(h(i));
                end % if i==10

            end % for i = 1:length(block_names)

            sgtitle(sprintf('%s %s %s %s %s',sub_name, task_name, metric_name, band_name, measure_name));

        else
            if (sum(isnan(tab2use.(task_name)))==length(tab2use.(task_name)))
                data2plot = repmat(0, length(tab2use.(task_name)),1);
            else
                data2plot = tab2use.(task_name);
            end
            colorlimits = prctile(data2plot,colorlimits_percentiles);
            topoplot(data2plot, chanlocs, 'colormap',colormap2use,'maplimits',colorlimits); colorbar;
            title(sprintf('%s %s %s  %s %s',sub_name, task_name, metric_name, band_name, measure_name));

        end % if strcmp(task_name,'RestingState')

        % save figure
        saveas(gcf, fname2save);

    end % for iband
end % for imeas

% aperiodic measures --------------------------------------------------
freq_bands = {'aperiodic'};
measure_names = {'offset','exponent'};
for iband = 1:length(freq_bands)
    band_name = freq_bands{iband};
    for imeas = 1:length(measure_names)
        measure_name = measure_names{imeas};
        tab2use = fooof_results.(sprintf('%s_%s',band_name, measure_name));

        out_path = fileparts(out_fig_file);
        fname2save = fullfile(out_path, sprintf('%s_%s_%s_%s_%s_topoplot.pdf', ...
            sub_name, task_name, metric_name, ...
            band_name,measure_name));

        figure; set(gcf,'color','white');
        colormap2use = flipud(lbmap(100,'RedBlue'));
        colorlimits_percentiles = [10,90]; % percentiles for color scaling

        if strcmp(task_name,'RestingState')

            block_names = {'open1','open2','open3','open4','open5','closed1','closed2','closed3','closed4','closed5'};
            concat_data = [];
            for i = 1:length(block_names)
                concat_data = [concat_data; tab2use.(block_names{i})];
            end % for i = 1:length(block_names)
            colorlimits = prctile(concat_data,colorlimits_percentiles);
            tlo = tiledlayout(2,length(block_names)/2);

            for i = 1:length(block_names)
                h(i) = nexttile(tlo);
                if (sum(isnan(tab2use.(block_names{i})))==length(tab2use.(block_names{i})))
                    data2plot = repmat(0, length(tab2use.(block_names{i})),1);
                else
                    data2plot = tab2use.(block_names{i});
                end
                topoplot(data2plot, chanlocs, 'colormap', colormap2use, 'maplimits', colorlimits);
                title(block_names{i});

                if i==length(block_names)
                    cbh = colorbar(h(i));
                end % if i==10

            end % for i = 1:length(block_names)

            sgtitle(sprintf('%s %s %s %s %s',sub_name, task_name, metric_name, band_name, measure_name));

        else
            if (sum(isnan(tab2use.(task_name)))==length(tab2use.(task_name)))
                data2plot = repmat(0, length(tab2use.(task_name)),1);
            else
                data2plot = tab2use.(task_name);
            end
            colorlimits = prctile(data2plot,colorlimits_percentiles);
            topoplot(data2plot, chanlocs, 'colormap',colormap2use,'maplimits',colorlimits); colorbar;
            title(sprintf('%s %s %s  %s %s',sub_name, task_name, metric_name, band_name, measure_name));

        end % if strcmp(task_name,'RestingState')

        % save figure
        saveas(gcf, fname2save);

    end % for iband
end % for imeas

end % function cmi_make_topoplot_fooof


%% function for making the topoplot
function cmi_make_topoplot_other(data2use, chanlocs, task_name, sub_name, metric_name, out_fig_file)

figure; set(gcf,'color','white');
colormap2use = flipud(lbmap(100,'RedBlue'));
colorlimits_percentiles = [10,90]; % percentiles for color scaling

if strcmp(task_name,'RestingState')

    block_names = {'open1','open2','open3','open4','open5','closed1','closed2','closed3','closed4','closed5'};
    concat_data = [];
    for i = 1:length(block_names)
        concat_data = [concat_data; data2use.(block_names{i})];
    end % for i = 1:length(block_names)
    colorlimits = prctile(concat_data,colorlimits_percentiles);
    tlo = tiledlayout(2,length(block_names)/2);

    for i = 1:length(block_names)
        h(i) = nexttile(tlo);
        if (sum(isnan(data2use.(block_names{i})))==length(data2use.(block_names{i})))
            data2plot = repmat(0, length(data2use.(block_names{i})),1);
        else
            data2plot = data2use.(block_names{i});
        end
        topoplot(data2plot, chanlocs, 'colormap', colormap2use, 'maplimits', colorlimits);
        title(block_names{i});

        if i==length(block_names)
            cbh = colorbar(h(i));
        end % if i==10

    end % for i = 1:length(block_names)

    sgtitle(sprintf('%s %s %s',sub_name, task_name, metric_name));

else

    if (sum(isnan(data2use.(task_name)))==length(data2use.(task_name)))
        data2plot = repmat(0, length(data2use.(task_name)),1);
    else
        data2plot = data2use.(task_name);
    end
    colorlimits = prctile(data2plot,colorlimits_percentiles);
    topoplot(data2plot, chanlocs, 'colormap',colormap2use,'maplimits',colorlimits); colorbar;
    title(sprintf('%s %s %s',sub_name, task_name, metric_name));

end % if strcmp(task_name,'RestingState')

% save figure
saveas(gcf, out_fig_file);

end % function cmi_make_topoplot_other




%
% %% function for making the topoplot
% function cmi_make_topoplot(data2use, chanlocs, task_name, sub_name, metric_name, out_fig_file)
%
% if strcmp(metric_name,'fooof')
%
%     fooof_results = cmi_reduce_fooof(data2use);
%
%     if strcmp(task_name,'RestingState')
%         block_names = {'open1','open2','open3','open4','open5','closed1','closed2','closed3','closed4','closed5'};
%     else
%         block_names = task_name;
%     end % if strcmp(task_name,'RestingState')
%
%     % oscillations in different frequency bands ---------------------------
%     freq_bands = {'delta','theta','alpha','beta','gamma'};
%     measure_names = {'freq','amp','bandwidth'};
%     for iband = 1:length(freq_bands)
%         band_name = freq_bands{iband};
%         for imeas = 1:length(measure_names)
%             measure_name = measure_names{imeas};
%             tab2use = fooof_results.(sprintf('%s_%s',band_name, measure_name));
%
%             out_path = fileparts(out_fig_file);
%             fname2save = fullfile(out_path, sprintf('%s_%s_%s_%s_%s_topoplot.pdf', ...
%                 sub_name, task_name, metric_name, ...
%                 band_name,measure_name));
%
%             %             cmi_fooof_plot(tab2use, block_names, chanlocs, ...
%             %                 sub_name, task_name, metric_name, ...
%             %                 band_name, measure_name, fname2save);
%
%             % -------------------------------------------------------------------------
%             figure; set(gcf,'color','white');
%             colormap2use = flipud(lbmap(100,'RedBlue'));
%             colorlimits_percentiles = [10,90]; % percentiles for color scaling
%
%             if strcmp(task_name,'RestingState')
%
%                 block_names = {'open1','open2','open3','open4','open5','closed1','closed2','closed3','closed4','closed5'};
%                 concat_data = [];
%                 for i = 1:length(block_names)
%                     concat_data = [concat_data; tab2use.(block_names{i})];
%                 end % for i = 1:length(block_names)
%                 colorlimits = prctile(concat_data,colorlimits_percentiles);
%                 tlo = tiledlayout(2,length(block_names)/2);
%
%                 for i = 1:length(block_names)
%                     h(i) = nexttile(tlo);
%                     topoplot(tab2use.(block_names{i}), chanlocs, 'colormap', colormap2use, 'maplimits', colorlimits);
%                     title(block_names{i});
%
%                     if i==length(block_names)
%                         cbh = colorbar(h(i));
%                     end % if i==10
%
%                 end % for i = 1:length(block_names)
%
%                 sgtitle(sprintf('%s %s %s %s %s',sub_name, task_name, metric_name, band_name, measure_name));
%
%             else
%
%                 data2plot = tab2use.(task_name);
%                 colorlimits = prctile(data2plot,colorlimits_percentiles);
%                 topoplot(data2plot, chanlocs, 'colormap',colormap2use,'maplimits',colorlimits); colorbar;
%                 title(sprintf('%s %s %s  %s %s',sub_name, task_name, metric_name, band_name, measure_name));
%
%             end % if strcmp(task_name,'RestingState')
%
%             % save figure
%             saveas(gcf, fname2save);
%             %--------------------------------------------------------------------------
%
%         end % for iband
%     end % for imeas
%
%     % aperiodic measures --------------------------------------------------
%     freq_bands = {'aperiodic'};
%     measure_names = {'offset','exponent'};
%     for iband = 1:length(freq_bands)
%         band_name = freq_bands{iband};
%         for imeas = 1:length(measure_names)
%             measure_name = measure_names{imeas};
%             tab2use = fooof_results.(sprintf('%s_%s',band_name, measure_name));
%
%             out_path = fileparts(out_fig_file);
%             fname2save = fullfile(out_path, sprintf('%s_%s_%s_%s_%s_topoplot.pdf', ...
%                 sub_name, task_name, metric_name, ...
%                 band_name,measure_name));
%
%             %             cmi_fooof_plot(tab2use, block_names, chanlocs, ...
%             %                 sub_name, task_name, metric_name, ...
%             %                 band_name, measure_name, fname2save);
%
%             % -------------------------------------------------------------------------
%             figure; set(gcf,'color','white');
%             colormap2use = flipud(lbmap(100,'RedBlue'));
%             colorlimits_percentiles = [10,90]; % percentiles for color scaling
%
%             if strcmp(task_name,'RestingState')
%
%                 block_names = {'open1','open2','open3','open4','open5','closed1','closed2','closed3','closed4','closed5'};
%                 concat_data = [];
%                 for i = 1:length(block_names)
%                     concat_data = [concat_data; tab2use.(block_names{i})];
%                 end % for i = 1:length(block_names)
%                 colorlimits = prctile(concat_data,colorlimits_percentiles);
%                 tlo = tiledlayout(2,length(block_names)/2);
%
%                 for i = 1:length(block_names)
%                     h(i) = nexttile(tlo);
%                     topoplot(tab2use.(block_names{i}), chanlocs, 'colormap', colormap2use, 'maplimits', colorlimits);
%                     title(block_names{i});
%
%                     if i==length(block_names)
%                         cbh = colorbar(h(i));
%                     end % if i==10
%
%                 end % for i = 1:length(block_names)
%
%                 sgtitle(sprintf('%s %s %s %s %s',sub_name, task_name, metric_name, band_name, measure_name));
%
%             else
%
%                 data2plot = tab2use.(task_name);
%                 colorlimits = prctile(data2plot,colorlimits_percentiles);
%                 topoplot(data2plot, chanlocs, 'colormap',colormap2use,'maplimits',colorlimits); colorbar;
%                 title(sprintf('%s %s %s  %s %s',sub_name, task_name, metric_name, band_name, measure_name));
%
%             end % if strcmp(task_name,'RestingState')
%
%             % save figure
%             saveas(gcf, fname2save);
%             %--------------------------------------------------------------------------
%
%         end % for iband
%     end % for imeas
%
% else
%
%     figure; set(gcf,'color','white');
%     colormap2use = flipud(lbmap(100,'RedBlue'));
%     colorlimits_percentiles = [10,90]; % percentiles for color scaling
%
%     if strcmp(task_name,'RestingState')
%
%         block_names = {'open1','open2','open3','open4','open5','closed1','closed2','closed3','closed4','closed5'};
%         concat_data = [];
%         for i = 1:length(block_names)
%             concat_data = [concat_data; data2use.(block_names{i})];
%         end % for i = 1:length(block_names)
%         colorlimits = prctile(concat_data,colorlimits_percentiles);
%         tlo = tiledlayout(2,length(block_names)/2);
%
%         for i = 1:length(block_names)
%             h(i) = nexttile(tlo);
%             topoplot(data2use.(block_names{i}), chanlocs, 'colormap', colormap2use, 'maplimits', colorlimits);
%             title(block_names{i});
%
%             if i==length(block_names)
%                 cbh = colorbar(h(i));
%             end % if i==10
%
%         end % for i = 1:length(block_names)
%
%         sgtitle(sprintf('%s %s %s',sub_name, task_name, metric_name));
%
%     else
%
%         data2plot = data2use.(task_name);
%         colorlimits = prctile(data2plot,colorlimits_percentiles);
%         topoplot(data2plot, chanlocs, 'colormap',colormap2use,'maplimits',colorlimits); colorbar;
%         title(sprintf('%s %s %s',sub_name, task_name, metric_name));
%
%     end % if strcmp(task_name,'RestingState')
%
%     % save figure
%     saveas(gcf, out_fig_file);
%
% end % if strcmp(metric_name,'fooof')
%
%
%
%
% end % function cmi_make_topoplot


%% function lbmap
function map = lbmap(n,scheme)

%defensive programming
error(nargchk(0,2,nargin))
error(nargoutchk(0,1,nargout))

%defaults
if nargin<2
  scheme = 'Blue';
end
if nargin<1
  n = size(get(gcf,'colormap'),1);
end

%valid schemes
switch lower(scheme)
  case 'blue'
    baseMap = BlueMap;
  case 'bluegray'
    baseMap = BlueGrayMap;
  case 'brownblue'
    baseMap = BrownBlueMap;
  case 'redblue'
    baseMap = RedBlueMap;
  otherwise
    error(['Invalid scheme ' scheme])
end
idx1 = linspace(0,1,size(baseMap,1));
idx2 = linspace(0,1,n);
map = interp1(idx1,baseMap,idx2);

end % function lbmap

function baseMap = BlueMap
baseMap = [243 246 248;
           224 232 240;
           171 209 236;
           115 180 224;
            35 157 213;
             0 142 205;
             0 122 192]/255;
end % function blueMap

function baseMap = BlueGrayMap
%DivergingBlueGray
baseMap = [  0 170 227;
            53 196 238;
           133 212 234;
           190 230 242;
           217 224 230;
           146 161 170;
           109 122 129;
            65  79  81]/255;
end % function BlueGrayMap

function baseMap = BrownBlueMap
baseMap = [144 100  44;
           187 120  54;
           225 146  65;
           248 184 139;
           244 218 200;
           241 244 245;
           207 226 240;
           160 190 225;
           109 153 206;
            70  99 174;
            24  79 162]/255;
end % function BrownBlueMap

function baseMap = RedBlueMap
baseMap = [175  53  71;
           216  82  88;
           239 133 122;
           245 177 139;
           249 216 168;
           242 238 197;
           216 236 241;
           154 217 238;
            68 199 239;
             0 170 226;
             0 116 188]/255;
end % function RedBlueMap
