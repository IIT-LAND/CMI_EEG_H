function cmi_tsplot(datafile)
% cmi_tsplot - make a plot of the timeseries
%
%   INPUT
%       datafile = full filename to the *.set file you want to plot
%
%   Example usage:
%       datafile = '/media/DATA/RAW/cmihbn/data/raw/NDARVB151GXF/Present/NDARVB151GXF_Present.set'
%       run_on_server = true;
%       cmi_tsplot(datafile, run_on_server);
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

try
    winlength = round((EEG.times(end)/1000)-2);
    eegplot('noui', EEG.data,'winlength',winlength, 'color','on', ...
        'title',plot_title, 'plottitle',plot_title);
    set(gcf,'color','white');
catch
    figure; set(gcf,'color','white');
    plot(EEG.data'); grid on;
    title(plot_title);
end

%% save the plot
outfile = fullfile(fpath,sprintf('%s_%s_%s_tsplot.jpg',sub_name,task_name,dtype));
% saveas(gcf,outfile);
print(gcf,outfile, '-djpeg','-noui');
close(gcf);

end % function cmi_tsplot


