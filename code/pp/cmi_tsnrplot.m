function cmi_tsnrplot(datafile)
% cmi_snrplot - make a temporal signal-to-noise (tSNR) plot 
%
%   INPUT
%       datafile = full filename to the *.set file you want to plot
%       outfile = full filename to the jpg plot you want to save
%
%   Example usage:
%       datafile = '/media/DATA/RAW/cmihbn/data/raw/NDARVB151GXF/Present/NDARVB151GXF_Present.set'
%       run_on_server = true;
%       cmi_tsnrplot(datafile, run_on_server);
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

%% compute tSNR
electrode_mean = nanmean(EEG.data,2);
t_mean = repmat(electrode_mean,1,size(EEG.data,2));
electrode_sd = nanstd(EEG.data,0,2);
t_sd = repmat(electrode_sd,1,size(EEG.data,2));
tsnr_mat = (EEG.data - t_mean)./t_sd;
% tsnr_electrode = electrode_mean./electrode_sd;
tsnr_electrode = (abs(electrode_mean./electrode_sd)).*100;
grand_mean_tsnr = nanmean(tsnr_electrode);

%% make tSNR carpet plot
plot_title = sprintf('%s %s tSNR',sub_name,task_name);

figure; set(gcf,'color','white');
imagesc(tsnr_mat); colormap(gray); colorbar;
xlabel('Time'); ylabel('Electrodes');
title(plot_title);

% save the plot
outfile = fullfile(fpath,sprintf('%s_%s_%s_tsnrcarpetplot.jpg',sub_name,task_name,dtype));
% saveas(gcf,outfile);
print(gcf,outfile, '-djpeg','-noui');
close(gcf);

%% make tSNR topoplot
plot_title = sprintf('%s %s tSNR',sub_name,task_name);
colormap2use = flipud(lbmap(100,'RedBlue'));

figure; set(gcf,'color','white');
colorlimits_percentiles = [10,90]; % percentiles for color scaling
colorlimits = prctile(tsnr_electrode,colorlimits_percentiles);
% colorlimits = [0,1];
topoplot(tsnr_electrode, EEG.chanlocs, 'colormap', colormap2use, 'maplimits', colorlimits); colorbar;
title(plot_title);

% save the plot
outfile = fullfile(fpath,sprintf('%s_%s_%s_tsnrtopoplot.jpg',sub_name,task_name,dtype));
% saveas(gcf,outfile);
print(gcf,outfile, '-djpeg','-noui');
close(gcf);


%% save out the mean tSNR over all electrodes
if strcmp(dtype,'preproc')
    tab2use = array2table(grand_mean_tsnr,'VariableNames',{'mean_tSNR'});
    writetable(tab2use, fullfile(fpath,sprintf('%s_%s_%s_mean_tsnr.csv',sub_name,task_name,dtype)));
end % if strcmp(dtype,'preproc')

end % function cmi_tsnrplot




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
