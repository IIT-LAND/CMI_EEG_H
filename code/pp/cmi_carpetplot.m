function cmi_carpetplot(datafile)
% cmi_carpetplot - make a carpet plot of the timeseries
%
%   INPUT
%       datafile = full filename to the *.set file you want to plot
%       outfile = full filename to the jpg plot you want to save
%
%   Example usage:
%       datafile = '/media/DATA/RAW/cmihbn/data/raw/NDARVB151GXF/Present/NDARVB151GXF_Present.set'
%       outfile = '/media/DATA/RAW/cmihbn/data/raw/NDARVB151GXF/Present/NDARVB151GXF_Present_raw_carpetplot.jpg'
%       run_on_server = true;
%       cmi_carpetplot(datafile, outfile, run_on_server);
%
%   written by mvlombardo
%

%% parse datafile
[fpath] = fileparts(datafile);
[subpath, task_name] = fileparts(fpath);
[dtypepath, sub_name] = fileparts(subpath);
[datapath, dtype] = fileparts(dtypepath);
[rootpath] = fileparts(datapath);
codepath = fullfile(rootpath,'code');

if strcmp(rootpath,'/media/DATA/RAW/cmihbn')
    eeglab_dir = '/home/mlombardo/eeglab_20201226';
else
    eeglab_dir = '/Users/mlombardo/Dropbox/matlab/eeglab_20201226';
end
cd(eeglab_dir);
eeglab nogui;
cd(codepath);

%% load data
EEG = pop_loadset(datafile);

%% make the plot
plot_title = sprintf('%s %s',sub_name,task_name);

figure; set(gcf,'color','white');
imagesc(EEG.data); colormap(gray); colorbar;
xlabel('Time'); ylabel('Electrodes');
title(plot_title);

%% save the plot
outfile = fullfile(fpath,sprintf('%s_%s_%s_carpetplot.jpg',sub_name,task_name,dtype));
% saveas(gcf,outfile);
print(gcf,outfile, '-djpeg','-noui');
close(gcf);

end % function cmi_carpetplot


