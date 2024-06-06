function results = cmi_preproc_quality(rawfile,ppfile)
% cmi_preproc_quality
%
%   Compute measures from Automagic toolbox that quantify quality metrics
%   for a dataset.
%
%   INPUT
%       
%       rawfile = a raw *.set file 
%       ppfile = a preproc *.set file
%

%% parameters
overallThresh = 20; 
timeThresh = 10;
chanThresh = 10;

%% Parse apart information from filenames
[fpath] = fileparts(ppfile);
[subpath,task_name] = fileparts(fpath);
[preprocpath,sub_name] = fileparts(subpath);
[datapath] = fileparts(preprocpath);
[rootpath] = fileparts(datapath);
codepath = fullfile(rootpath,'code');

if strcmp(rootpath,'/media/DATA/RAW/cmihbn')
    eeglab_dir = '/home/mlombardo/eeglab_20201226';
    automagic_dir = '/home/mlombardo/automagic';
else
    eeglab_dir = '/Users/mlombardo/Dropbox/matlab/eeglab_20201226';
    automagic_dir = '/Users/mlombardo/Dropbox/matlab/automagic';
end
cd(eeglab_dir);
eeglab nogui;
cd(codepath);
cd(automagic_dir);
addAutomagicPaths;
cd(codepath);

%% load data
EEGraw = pop_loadset(rawfile);
EEGpp = pop_loadset(ppfile);

%% calculate quality metrics
try
    Qraw = calcQuality(EEGraw,[], ...
        'overallThresh',overallThresh, ...
        'timeThresh',timeThresh, ...
        'chanThresh',chanThresh);
catch
    Qraw.OHA = NaN;
    Qraw.THV = NaN;
    Qraw.CHV = NaN;
    Qraw.MAV = NaN;
    Qraw.RBC = NaN;
    Qraw.settings.avReg = 1;
    Qraw.settings.chanThresh = chanThresh;
    Qraw.settings.checkboxCutoff_CHV = 0;
    Qraw.settings.Cutoff_CHV = 100;
    Qraw.settings.overallThresh = overallThresh;
    Qraw.settings.RejRatio_CHV = 0.5;
    Qraw.settings.timeThresh = timeThresh;
end % try

try
    Qpp = calcQuality(EEGpp,[], ...
        'overallThresh',overallThresh, ...
        'timeThresh',timeThresh, ...
        'chanThresh',chanThresh);
catch
    Qpp.OHA = NaN;
    Qpp.THV = NaN;
    Qpp.CHV = NaN;
    Qpp.MAV = NaN;
    Qpp.RBC = NaN;
    Qpp.settings.avReg = 1;
    Qpp.settings.chanThresh = chanThresh;
    Qpp.settings.checkboxCutoff_CHV = 0;
    Qpp.settings.Cutoff_CHV = 100;
    Qpp.settings.overallThresh = overallThresh;
    Qpp.settings.RejRatio_CHV = 0.5;
    Qpp.settings.timeThresh = timeThresh;
end

%% save metrics
data2use = [Qraw.OHA,Qpp.OHA,Qraw.THV,Qpp.THV,Qraw.CHV,Qpp.CHV,Qraw.MAV,Qpp.MAV];
tab2use = array2table(data2use,'VariableNames',{'raw_OHA','preproc_OHA', ...
    'raw_THV','preproc_THV', ...
    'raw_CHV','preproc_CHV', ...
    'raw_MAV','preproc_MAV'});

fname2save = fullfile(fpath,sprintf('%s_%s_preproc_dataQuality.csv',sub_name,task_name));
writetable(tab2use,fname2save);

%% plot metrics
figure; set(gcf,'color','white');

% OHV
subplot(2,2,1);
bar([tab2use.raw_OHA,tab2use.preproc_OHA],'FaceColor',[0.5,0.5,0.5]);
hold on;
plot([tab2use.raw_OHA,tab2use.preproc_OHA],'k','LineWidth',2);
scatter(1:2,[tab2use.raw_OHA,tab2use.preproc_OHA],'k');
ylim([-0.05,1.05]); grid on;
yLims = ylim;
text(1.75,(yLims(2)-yLims(2)*0.1),sprintf('OHA = %0.4f',tab2use.preproc_OHA));
ylabel('Overall High Amplitude (OHA)');
xticklabels({'Raw','Preproc'});
title(sprintf('Ratio of data points > %d mV (OHV)',overallThresh))

% THV
subplot(2,2,2);
bar([tab2use.raw_THV,tab2use.preproc_THV],'FaceColor',[0.5,0.5,0.5]);
hold on;
plot([tab2use.raw_THV,tab2use.preproc_THV],'k','LineWidth',2);
scatter(1:2,[tab2use.raw_THV,tab2use.preproc_THV],'k');
ylim([-0.05,1.05]); grid on;
yLims = ylim;
text(1.75,(yLims(2)-yLims(2)*0.1),sprintf('THV = %0.4f',tab2use.preproc_THV));
ylabel('Timepoints of High Variance (THV)');
xticklabels({'Raw','Preproc'});
title(sprintf('Ratio of time points SD > %d mV (THV)',timeThresh))

% CHV
subplot(2,2,3);
bar([tab2use.raw_CHV,tab2use.preproc_CHV],'FaceColor',[0.5,0.5,0.5]);
hold on;
plot([tab2use.raw_CHV,tab2use.preproc_CHV],'k','LineWidth',2);
scatter(1:2,[tab2use.raw_CHV,tab2use.preproc_CHV],'k');
ylim([-0.05,1.05]); grid on;
yLims = ylim;
text(1.75,(yLims(2)-yLims(2)*0.1),sprintf('CHV = %0.4f',tab2use.preproc_CHV));
ylabel('Channels of High Variance (CHV)');
xticklabels({'Raw','Preproc'});
title(sprintf('Ratio of channels SD > %d mV (THV)',chanThresh))

% MAV
subplot(2,2,4);
bar([tab2use.raw_MAV,tab2use.preproc_MAV],'FaceColor',[0.5,0.5,0.5]);
hold on;
plot([tab2use.raw_MAV,tab2use.preproc_MAV],'k','LineWidth',2);
scatter(1:2,[tab2use.raw_MAV,tab2use.preproc_MAV],'k');
% ylim([-0.05,1.05]);
grid on;
yLims = ylim;
text(1.75,(yLims(2)-yLims(2)*0.1),sprintf('MAV = %0.4f',tab2use.preproc_MAV));
ylabel('Mean Absolute Voltage (MAV)');
xticklabels({'Raw','Preproc'});
title('Mean Absolute Voltage (MAV)')

sgtitle(sprintf('%s %s Data Quality Metrics',sub_name,task_name));

fname2save = fullfile(fpath,sprintf('%s_%s_preproc_dataQuality.jpg',sub_name,task_name));
print(gcf,fname2save, '-djpeg','-noui');


% end function cmi_preproc_quality