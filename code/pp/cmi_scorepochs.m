function [idx_best_ep, epoch, score_Xep, score_chXep, psd_chXep, pxx] = cmi_scorepochs(datafile, cfg) 
% cmi_scorepochs 
%
% Function to select the best (most homogenoous) M/EEG epochs from a
% resting-state recordings. 
%
%     Copyright (C) 2020 Matteo Demuru, Matteo Fraschini
%     last modified by andrea.vitale@gmail.com  20210706
%
% INPUT
%    datafile            - file .set (as in the eeglab structure with already the channels location) 
%    cfg struct with the following fields
%           freqRange    - array with the frequency range used to compute the power
%                          spectrum (see MATLAB pwelch function)
%           fs           - integer representing sample frequency         
%           windowL      - integer representing the window length (in seconds)  
%           smoothFactor - smoothing factor for the power spectrum
%     
% OUTPUT
%      
%    epoch       -  cell array of the data divided in equal length epochs 
%                   of length windowL (channels X time samples)
%                  
%    idx_best_ep - array of indexes sorted according to the best score
%                  this array should be used for the selection of the best
%                  epochs
%
%    score_Xep   - array of score per epoch
%
%    ADDITIONAL Output:
%           score_chXep:  scorepoch values not averaged across channels
%
%           psd_chXep: psd computed for each epoch and channel
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <https://www.gnu.org/licenses/>.


%% parse input file path and names
[fpath] = fileparts(datafile);
[subpath, task_name] = fileparts(fpath);
[preprocpath, sub_name] = fileparts(subpath);
[datapath] = fileparts(preprocpath);
[rootpath] = fileparts(datapath);
codepath = fullfile(rootpath,'code');

%% Install eeglab
if strcmp(rootpath,'/media/DATA/RAW/cmihbn')
    eeglab_dir = '/home/mlombardo/eeglab_20201226';
else
    eeglab_dir = '/Users/mlombardo/Dropbox/matlab/eeglab_20201226';
end % if strcmp(rootpath,'/media/DATA/RAW/cmihbn')
cd(eeglab_dir);
eeglab nogui;
cd(codepath);

%% load dataset
try
    eeg_struct = pop_loadset(datafile);
catch ME
    disp(ME.message);
end


%% configuration 
if ~exist('cfg') || isempty('cfg')
    
    cfg = [];
    
    % frequency range
    cfg.freqRange = [2:50];
    
    % sampling rate
    cfg.fs = eeg_struct.srate;
    
    % window length in seconds
    cfg.windowL = 2; 
    
    % smoothing factor
    cfg.smoothFactor = 0;
    
    % do plotting
    cfg.do_plot = true;
    
    % save results to a file
    cfg.save_results = true;
    
    % file stems to put on the figures
    cfg.fig1_fname_fstem = 'scorepochs1';
    cfg.fig2_fname_fstem = 'scorepochs2';
    
end % if ~exist('cfg') || isempty('cfg')


%% divide the data in epochs of windowL length (in sec)
data = eeg_struct.data;

epLen = cfg.windowL * cfg.fs;           % epoch length in samples
dataLen = size(data,2);                 % number of samples
nCh = size(data,1);                     % number of channels
idx_ep = 1:epLen:(dataLen - epLen + 1); % start index for epochs
nEp = numel(idx_ep);                    % number of epochs

% pre-allocate cell arrays
epoch   = cell(1,nEp);
pxx     = cell(1,nEp); % for storing the results of pwelch(epoch_i)

% loop over epochs
for i_epoch = 1:nEp
    
    % grab data for the specific epoch of interest
    samples2use = idx_ep(i_epoch):(idx_ep(i_epoch)+epLen-1);
    epoch{i_epoch} = data(:,samples2use);
    
    % compute POWER SPECTRUM
    % pxx{e} dimensions = n_channel x n_freq
    % (power for each frequency in a specific epoch)
    pxx{i_epoch} = pwelch(epoch{i_epoch}',[],[], cfg.freqRange, cfg.fs)';
    
    % apply smoothing factor if necessary
    if(cfg.smoothFactor ~= 0)
        pxx{i_epoch} = movmean(pxx{i_epoch}',cfg.smoothFactor)';
    end % if(cfg.smoothFactor ~= 0)
    
end % for i_epoch = 1:nEp

%% compute SCORE across channels and across epochs

score_chXep = zeros(nCh,nEp); % pre-allocate score array [channels, epochs]
psd_chXep = zeros(nCh, nEp, numel(cfg.freqRange)); % pre-allocate empty array

% loop over channels
for i_channel = 1:nCh
    
    % loop over epochs
    for i_epoch = 1:nEp
        
        % grab PSDs as an array of [channels, epochs, frequencies]
        psd_chXep(i_channel,i_epoch,:) =  pxx{i_epoch}(i_channel,:);
        
    end % for i_epoch = 1:nEp
    
    % compute scores as spearman correlation
    data2use = squeeze(psd_chXep(i_channel,:,:))';
    % [n_epochs, n_epochs] correlation matrix for the channel of interest
    score_ch = corr(data2use,'type','Spearman'); 
    
    % computes the mean of each column in the [n_epochs,n_epochs] correlation matrix
    score_chXep(i_channel,:) = mean(score_ch); 
    
end % for i_channel = 1:nCh

% mean over all the channels, resulting in one score per each epoch, 
% representing the mean score for that epoch over all channels
score_Xep = mean(score_chXep,1); 

[~,idx_best_ep] = sort(score_Xep,'descend');


%% Plots
if cfg.do_plot
    
    % FIGURE 1: BAR GRAPH of the score_epoch values (not sorted)
    thresh_level = 2; %std deviation
    line_width = 1.5;
    
    figure; set(gcf,'color','white');
    subplot(2,2,1); 
    hold on; 
    pl1 = bar(score_Xep);
    set(pl1,'FaceColor',[ 1 1 1 ]);
    xlim([1 length(score_Xep)]);
    xlabel('Epoch');
    ylabel('Score');
    ylim([0 1.05]);
    grid on;    
    yline(mean(score_Xep),'k-', 'LineWidth',line_width);
    yline(mean(score_Xep) - std(score_Xep),'k--', 'LineWidth',line_width);
    yline(mean(score_Xep) - 2*std(score_Xep),'k:', 'LineWidth',line_width);
    title(sprintf('Average Epoch Score = %.2f',mean(score_Xep)));
    
    % identify CHANNELS with low score_epoch as compared to other channels:
    subplot(2,2,2); 
    hold on;
    
    score_xchan_mean = mean(score_chXep,2);
    score_xchan_std = std(score_chXep,1,2);
    score_allchan_mean = mean(score_xchan_mean,1);
    score_allchan_std = std(score_xchan_mean,1,1);
    errorbar(score_xchan_mean, score_xchan_std,'ko');

    hold on;
    yline(score_allchan_mean, 'k','LineWidth',line_width);
    yline(score_allchan_mean - score_allchan_std,'k--','LineWidth',line_width);
    yline(score_allchan_mean - 2*score_allchan_std, 'k:','LineWidth',line_width);
    xlabel('Channels');
    ylabel('Score');
%     ylim([0 1.05]);
    grid on;
    lower_bound_y = (min(score_Xep(:)) - max(score_xchan_std(:))) - 0.1;
    ylim([lower_bound_y, 1.05]);
    
    % find channel outliers:
    bad_chan_idx = find(score_xchan_mean <= (score_allchan_mean - thresh_level*score_allchan_std))';
    for ii = 1:length(bad_chan_idx)
        hold on;
        errorbar(bad_chan_idx(ii), score_xchan_mean(bad_chan_idx(ii)), ...
            score_xchan_std(bad_chan_idx(ii)),'rx', 'LineWidth', 2);
    end % for ii = 1:length(bad_chan_idx)
    
    title(sprintf('Bad Channels (Score <%d SD) in Red',thresh_level));
%     title(['Bad Channels (Score < ' num2str(thresh_level) 'SD) = '...
%         num2str(reshape(bad_chan_idx,[1,length(bad_chan_idx)]))]);
    sgtitle(sprintf('%s %s',sub_name, task_name));
    
    % epoch with the LOWEST score
    [tmp, i_epoch] = min(score_Xep);  
    i_epoch = i_epoch(1);
    
    % change color of the histogram
    subplot(2,2,1); 
    hold on;
    bar(i_epoch,score_Xep(i_epoch),'r');
    ylim([lower_bound_y, 1.05]);

    subplot(2,2,3);
    imagesc(score_chXep); 
    colormap(flipud(lbmap(100,'RedBlue'))); colorbar;
    ylabel('Channels'); 
    xlabel('Epochs');
    title('Epoch Scores Per Channel');
    
    subplot(2,2,4);
    histogram(score_chXep(:),100, 'FaceColor',[0.5,0.5,0.5]); 
    grid on;
    ylabel('Count');
    xlabel('Score');
    title('All Epochs and Channels');    

    fname2save = fullfile(fpath,sprintf('%s_%s_preproc_%s.jpg', ...
        sub_name,task_name,cfg.fig1_fname_fstem));
    print(gcf,fname2save, '-djpeg','-noui');
    
    % MULTICHANNEL plot + TOPOPLOT
    % FIGURE 2 
    epoch_data = epoch{1,i_epoch};
    epoch_psd = pwelch(epoch_data,[],[],cfg.freqRange(1):cfg.freqRange(end),cfg.fs)';
    epoch_psd_xch = squeeze(psd_chXep(:,i_epoch,:));
    
    % threshold level for selecting NOISY CHANNELS:
    % +/- n  standard deviation
    thresh_level = 2;
    
    lowfreq_thresh = 10;
    highfreq_thresh = 40;
    freq_thresh = [ lowfreq_thresh  highfreq_thresh];
    
    chanloc = eeg_struct.chanlocs;
    n_chan = length(chanloc);
    
    % compare the psd of of all the other channels in a single epoch)
    % compute the psd mean chan by chan
    n_freq = size(epoch_psd,2);
    psd_mean = mean(epoch_psd_xch,1);
    psd_std = std(epoch_psd_xch,1);
    
    thresh_xfreq = thresh_level * psd_std; % noise threshold = +/- standard deviation
    
    % identify bad channels
    bad_chan_table = zeros(n_chan, n_freq);
    for i_chan = 1:n_chan

        for i_freq = 1:n_freq
        
            if epoch_psd_xch(i_chan, i_freq) > (psd_mean(i_freq)+thresh_xfreq(i_freq)) || ...
                    epoch_psd_xch(i_chan, i_freq) < (psd_mean(i_freq)-thresh_xfreq(i_freq))
                bad_chan_table(i_chan, i_freq) = 1;
            end % if epoch_psd_xch
            
        end % for i_freq = 1:n_freq
        
    end % for i_chan = 1:n_chan
    
    scroll_topoplot(epoch_data, bad_chan_table, chanloc, freq_thresh);
    title([' epoch ' num2str(i_epoch) ]);
    set(gcf,'color','white');
    
    fname2save = fullfile(fpath,sprintf('%s_%s_preproc_%s.jpg', ...
        sub_name,task_name,cfg.fig2_fname_fstem));
    print(gcf,fname2save, '-djpeg','-noui');
    
end % if cfg.do_plot

%% save out results to a file
if cfg.save_results

    for i = 1:length(score_Xep)
        epoch_names{i,1} = sprintf('epoch_%04d',i);
    end
    
    tab2write = cell2table([epoch_names num2cell(score_Xep)'],'VariableNames',{'Epoch','Score'});
    fname2save = fullfile(fpath,sprintf('%s_%s_preproc_scorepochs_avgEpochScores.csv',sub_name,task_name));
    writetable(tab2write,fname2save);
    
    tab2write = cell2table(num2cell(mean(score_Xep)),'VariableNames',{'Score'});
    fname2save = fullfile(fpath,sprintf('%s_%s_preproc_scorepochs_summaryEpochScore.csv',sub_name,task_name));
    writetable(tab2write,fname2save);
    
    chan_names = {chanloc.labels}';
    tab2write = cell2table([chan_names num2cell(score_chXep)],'VariableNames',[{'Channel'},epoch_names']);
    fname2save = fullfile(fpath,sprintf('%s_%s_preproc_scorepochs_chanEpochScores.csv',sub_name,task_name));
    writetable(tab2write,fname2save);

end % if cfg.save_results


fprintf('Done \n');

end % function cmi_scorepochs
    

%% function scroll_topoplot
function scroll_topoplot(epoch_data, bad_chan_table, chanloc, freq_thresh)
    %INPUT example:
    % chanloc = eeg_struct.chanlocs;
    
    lowfreq_thresh = freq_thresh(1);
    highfreq_thresh = freq_thresh(2);
    
    figure; 
    hold on;
    
    subplot(3,3,3); 
    hold on   
    marker_size = 20;
    
    for i_chan = 1:size(bad_chan_table,1)
        
        for i_freq = 1:size(bad_chan_table,2)
            
            if bad_chan_table(i_chan, i_freq) > 0
                
                if i_freq <= lowfreq_thresh
                    scatter(i_freq, i_chan, marker_size, 'k.');
                elseif i_freq >= highfreq_thresh
                    scatter(i_freq, i_chan, marker_size, 'r.');
                else
                    scatter(i_freq, i_chan, marker_size/5, 'b.');
                end % if i_freq <= lowfreq_thresh; scatter(i_freq, i_chan, marker_size, 'k.')
                
            end % if bad_chan_table(i_chan, i_freq) > 0
            
        end % for i_freq = 1:size(bad_chan_table,2)
        
    end % for i_chan = 1:size(bad_chan_table,1)
    xlim([ 1 size(bad_chan_table,2) ]) 
    ylim([ 1 size(bad_chan_table,1) ]) 
    ylabel('channels')
    xlabel('frequency')
        
    % LOW frequency (< 10Hz)
    bad_lowfreq_xchan = sum(bad_chan_table(:,1:lowfreq_thresh),2);
    subplot(3,3,6); hold on    
    %single vector of channel values
    topoplot(bad_lowfreq_xchan, chanloc, 'electrodes','on');
    colormap(flipud(gray));
    caxis([0 max(bad_lowfreq_xchan)]); %colorbar
    freezeColors;
    title(['chan outlier: LOW freq <' num2str(lowfreq_thresh) 'Hz'])
    
    % HIGH frequency (> 40Hz)
    bad_highfreq_xchan = sum(bad_chan_table(:,highfreq_thresh:end),2);
    subplot(3,3,9); hold on    
    topoplot(bad_highfreq_xchan, chanloc, 'electrodes','on');
    colormap(flipud(hot));
    caxis([0 max(bad_lowfreq_xchan)*2.5]); %colorbar
    freezeColors;
    title(['chan outlier: HIGH freq <' num2str(highfreq_thresh) 'Hz'])
    
    % MULTICHANNEL SCROLL 
    subplot(3,3,[1,2,4,5,7,8])
    plot_multichan_nonormalize(epoch_data, bad_chan_table); %colorbar
    %plot_multichan_nonormalize(epoch_data, bad_chan_table, freq_thresh); %colorbar
    xlabel('sample points')
    ylabel('channels')
    
end % function scroll_topoplot


%% function plot_multichan_nonormalize
function plot_multichan_nonormalize( x, bad_chan_table, freq_thresh )
%function plot_multichan( x, y, interval, normalize, bad_chan_idx )
%function plot_multichan( x, bad_chan_idx )

% Simple MATLAB function for plot multi-channel time-series data
% 
% Usage:
%    plot_multichan( y )    % <- y: signal
%    plot_multichan( x, y ) % <- x: time
% 
% Example: 
%    y = randn([20, 2000]); 
%    plot_multichan(y);
%
% Written by Hio-Been han, hiobeen.han@kaist.ac.kr, 2020-03-07
% modified by andrea.vitale@gmail.com  20210325

% parse arguments
if nargin==1, y=x; x = 1:length(y); end
y=x; x = 1:length(y); 
nChan = size(y,1);
if nChan > size(y,2),  y = y'; nChan = size(y,1); end
%if nargin < 4, normalize = 1; end
%if nargin < 4, normalize = 0; end
normalize = 0;
if normalize
    stds = nanstd( y, 0, 2 );
    for chIdx = 1:size(y,1), y(chIdx,:) = nanmean(stds) * (y(chIdx,:) / stds(chIdx)); end
end
% if nargin < 3
    interval = nanmean(range(y, 2)) * nChan / 2.5;
% end
y_center = linspace( -interval, interval, nChan );
y_center = linspace( -interval*1.15, interval*1.15, nChan );

if nargin < 3
    lowfreq_lim = 10; highfreq_lim = 65;
else
    lowfreq_lim = freq_thresh(1); 
    highfreq_lim = freq_thresh(2);
end

[bad_chan_idx, ~] = find(sum(bad_chan_table,2) >2 );
[bad_chan_lowfreq_idx, ~] = find(sum(bad_chan_table(:,1:lowfreq_lim),2) > 2);
[bad_chan_highfreq_idx, ~] = find(sum(bad_chan_table(:,highfreq_lim:end),2) > 2);

% set colormap
color_template =...
   [843 088 153;
    992 750 280;
    400 200 030;
    573 716 350;
    055 538 083]*.001;
c_space = repmat( color_template, [ ceil( nChan/size(color_template,1)), 1]);
% c_space=imresize(colormap('lines'),[nChan,3],'nearest');

% main plot
chanlab = {}; chanlab_pos = [];
lw = 0.5; %1;
chan_scale = 50; %micro_volt
for chanIdx = 1:nChan
    shift = y_center(chanIdx) + nanmean( y( chanIdx, : ), 2);
    
    if ismember(chanIdx, bad_chan_highfreq_idx) 
        plot( x, y( chanIdx, : ) - shift, 'r', 'LineWidth', lw*2); 
    %elseif ismember(chanIdx, bad_chan_idx) 
    elseif ismember(chanIdx, bad_chan_lowfreq_idx) 
        plot( x, y( chanIdx, : ) - shift, 'k', 'LineWidth', lw*2); 
    else
        plot( x, y( chanIdx, : ) - shift, 'b', 'LineWidth', lw/2);
        %plot( x, y( chanIdx, : ) - shift, 'Color', c_space( chanIdx,: ) , 'LineWidth', lw);
    end
    %ylim([ 0 chan_scale ])
    
%     if ismember(chanIdx, bad_chan_idx) 
%         chanIdx_reverse = nChan-chanIdx+1;
%         chanlab{chanIdx} = sprintf( 'Ch %02d', chanIdx_reverse);
%         chanlab_pos(chanIdx) =  y_center(chanIdx) ;
%     else
        chanIdx_reverse = nChan-chanIdx+1;
        if ismember(chanIdx_reverse, bad_chan_idx)
            chanlab{chanIdx} = sprintf('Ch %02d', chanIdx_reverse);
        else
            chanlab{chanIdx} = sprintf(' ');
        end
        chanlab_pos(chanIdx) =  y_center(chanIdx) ;
%    end
    
    
    if chanIdx ==1, hold on; end
end
hold off;

% enhance visibility
set(gca, 'YTick', chanlab_pos, 'YTickLabel', chanlab, 'FontSize', 5, ...
    'Clipping', 'on', 'Box', 'off', 'LineWidth', 2);
ylim([-1 1]*interval*1.2);
end % function plot_multichan_nonormalize
    

%% function freezeColors
function freezeColors(varargin)
% freezeColors  Lock colors of plot, enabling multiple colormaps per figure. (v2.3)
%
%   Problem: There is only one colormap per figure. This function provides
%       an easy solution when plots using different colomaps are desired 
%       in the same figure.
%
%   freezeColors freezes the colors of graphics objects in the current axis so 
%       that subsequent changes to the colormap (or caxis) will not change the
%       colors of these objects. freezeColors works on any graphics object 
%       with CData in indexed-color mode: surfaces, images, scattergroups, 
%       bargroups, patches, etc. It works by converting CData to true-color rgb
%       based on the colormap active at the time freezeColors is called.
%
%   The original indexed color data is saved, and can be restored using
%       unfreezeColors, making the plot once again subject to the colormap and
%       caxis.
%
%
%   Usage:
%       freezeColors        applies to all objects in current axis (gca),
%       freezeColors(axh)   same, but works on axis axh.
%
%   Example:
%       subplot(2,1,1); imagesc(X); colormap hot; freezeColors
%       subplot(2,1,2); imagesc(Y); colormap hsv; freezeColors etc...
%
%       Note: colorbars must also be frozen. Due to Matlab 'improvements' this can
%				no longer be done with freezeColors. Instead, please
%				use the function CBFREEZE by Carlos Adrian Vargas Aguilera
%				that can be downloaded from the MATLAB File Exchange
%				(http://www.mathworks.com/matlabcentral/fileexchange/24371)
%
%       h=colorbar; cbfreeze(h), or simply cbfreeze(colorbar)
%
%       For additional examples, see test/test_main.m
%
%   Side effect on render mode: freezeColors does not work with the painters
%       renderer, because Matlab doesn't support rgb color data in
%       painters mode. If the current renderer is painters, freezeColors
%       changes it to zbuffer. This may have unexpected effects on other aspects
%	      of your plots.
%
%       See also unfreezeColors, freezeColors_pub.html, cbfreeze.
%
%
%   John Iversen (iversen@nsi.edu) 3/23/05
%

%   Changes:
%   JRI (iversen@nsi.edu) 4/19/06   Correctly handles scaled integer cdata
%   JRI 9/1/06   should now handle all objects with cdata: images, surfaces, 
%                scatterplots. (v 2.1)
%   JRI 11/11/06 Preserves NaN colors. Hidden option (v 2.2, not uploaded)
%   JRI 3/17/07  Preserve caxis after freezing--maintains colorbar scale (v 2.3)
%   JRI 4/12/07  Check for painters mode as Matlab doesn't support rgb in it.
%   JRI 4/9/08   Fix preserving caxis for objects within hggroups (e.g. contourf)
%   JRI 4/7/10   Change documentation for colorbars

% Hidden option for NaN colors:
%   Missing data are often represented by NaN in the indexed color
%   data, which renders transparently. This transparency will be preserved
%   when freezing colors. If instead you wish such gaps to be filled with 
%   a real color, add 'nancolor',[r g b] to the end of the arguments. E.g. 
%   freezeColors('nancolor',[r g b]) or freezeColors(axh,'nancolor',[r g b]),
%   where [r g b] is a color vector. This works on images & pcolor, but not on
%   surfaces.
%   Thanks to Fabiano Busdraghi and Jody Klymak for the suggestions. Bugfixes 
%   attributed in the code.

% Free for all uses, but please retain the following:
%   Original Author:
%   John Iversen, 2005-10
%   john_iversen@post.harvard.edu

appdatacode = 'JRI__freezeColorsData';

[h, nancolor] = checkArgs(varargin);

%gather all children with scaled or indexed CData
cdatah = getCDataHandles(h);

%current colormap
cmap = colormap;
nColors = size(cmap,1);
cax = caxis;

% convert object color indexes into colormap to true-color data using 
%  current colormap
for hh = cdatah',
    g = get(hh);
    
    %preserve parent axis clim
    parentAx = getParentAxes(hh);
    originalClim = get(parentAx, 'clim');    
   
    %   Note: Special handling of patches: For some reason, setting
    %   cdata on patches created by bar() yields an error,
    %   so instead we'll set facevertexcdata instead for patches.
    if ~strcmp(g.Type,'patch'),
        cdata = g.CData;
    else
        cdata = g.FaceVertexCData; 
    end
    
    %get cdata mapping (most objects (except scattergroup) have it)
    if isfield(g,'CDataMapping'),
        scalemode = g.CDataMapping;
    else
        scalemode = 'scaled';
    end
    
    %save original indexed data for use with unfreezeColors
    siz = size(cdata);
    setappdata(hh, appdatacode, {cdata scalemode});

    %convert cdata to indexes into colormap
    if strcmp(scalemode,'scaled'),
        %4/19/06 JRI, Accommodate scaled display of integer cdata:
        %       in MATLAB, uint * double = uint, so must coerce cdata to double
        %       Thanks to O Yamashita for pointing this need out
        idx = ceil( (double(cdata) - cax(1)) / (cax(2)-cax(1)) * nColors);
    else %direct mapping
        idx = cdata;
        %10/8/09 in case direct data is non-int (e.g. image;freezeColors)
        % (Floor mimics how matlab converts data into colormap index.)
        % Thanks to D Armyr for the catch
        idx = floor(idx);
    end
    
    %clamp to [1, nColors]
    idx(idx<1) = 1;
    idx(idx>nColors) = nColors;

    %handle nans in idx
    nanmask = isnan(idx);
    idx(nanmask)=1; %temporarily replace w/ a valid colormap index

    %make true-color data--using current colormap
    realcolor = zeros(siz);
    for i = 1:3,
        c = cmap(idx,i);
        c = reshape(c,siz);
        c(nanmask) = nancolor(i); %restore Nan (or nancolor if specified)
        realcolor(:,:,i) = c;
    end
    
    %apply new true-color color data
    
    %true-color is not supported in painters renderer, so switch out of that
    if strcmp(get(gcf,'renderer'), 'painters'),
        set(gcf,'renderer','zbuffer');
    end
    
    %replace original CData with true-color data
    if ~strcmp(g.Type,'patch'),
        set(hh,'CData',realcolor);
    else
        set(hh,'faceVertexCData',permute(realcolor,[1 3 2]))
    end
    
    %restore clim (so colorbar will show correct limits)
    if ~isempty(parentAx),
        set(parentAx,'clim',originalClim)
    end
    
end %loop on indexed-color objects

end % function freezeColors


%% function getCDataHandles -- get handles of all descendents with indexed CData
function hout = getCDataHandles(h)
% getCDataHandles  Find all objects with indexed CData

%recursively descend object tree, finding objects with indexed CData
% An exception: don't include children of objects that themselves have CData:
%   for example, scattergroups are non-standard hggroups, with CData. Changing
%   such a group's CData automatically changes the CData of its children, 
%   (as well as the children's handles), so there's no need to act on them.

error(nargchk(1,1,nargin,'struct'))

hout = [];
if isempty(h),return;end

ch = get(h,'children');
for hh = ch'
    g = get(hh);
    if isfield(g,'CData'),     %does object have CData?
        %is it indexed/scaled?
        if ~isempty(g.CData) && isnumeric(g.CData) && size(g.CData,3)==1, 
            hout = [hout; hh]; %#ok<AGROW> %yes, add to list
        end
    else %no CData, see if object has any interesting children
            hout = [hout; getCDataHandles(hh)]; %#ok<AGROW>
    end
end

end % function getCDataHandles


%% function getParentAxes -- return handle of axes object to which a given object belongs
function hAx = getParentAxes(h)
% getParentAxes  Return enclosing axes of a given object (could be self)

error(nargchk(1,1,nargin,'struct'))
%object itself may be an axis
if strcmp(get(h,'type'),'axes'),
    hAx = h;
    return
end

parent = get(h,'parent');
if (strcmp(get(parent,'type'), 'axes')),
    hAx = parent;
else
    hAx = getParentAxes(parent);
end

end % function getParentAxes

%% function checkArgs -- Validate input arguments
function [h, nancolor] = checkArgs(args)
% checkArgs  Validate input arguments to freezeColors

nargs = length(args);
error(nargchk(0,3,nargs,'struct'))

%grab handle from first argument if we have an odd number of arguments
if mod(nargs,2),
    h = args{1};
    if ~ishandle(h),
        error('JRI:freezeColors:checkArgs:invalidHandle',...
            'The first argument must be a valid graphics handle (to an axis)')
    end
    % 4/2010 check if object to be frozen is a colorbar
    if strcmp(get(h,'Tag'),'Colorbar'),
      if ~exist('cbfreeze.m'),
        warning('JRI:freezeColors:checkArgs:cannotFreezeColorbar',...
            ['You seem to be attempting to freeze a colorbar. This no longer'...
            'works. Please read the help for freezeColors for the solution.'])
      else
        cbfreeze(h);
        return
      end
    end
    args{1} = [];
    nargs = nargs-1;
else
    h = gca;
end

%set nancolor if that option was specified
nancolor = [nan nan nan];
if nargs == 2,
    if strcmpi(args{end-1},'nancolor'),
        nancolor = args{end};
        if ~all(size(nancolor)==[1 3]),
            error('JRI:freezeColors:checkArgs:badColorArgument',...
                'nancolor must be [r g b] vector');
        end
        nancolor(nancolor>1) = 1; nancolor(nancolor<0) = 0;
    else
        error('JRI:freezeColors:checkArgs:unrecognizedOption',...
            'Unrecognized option (%s). Only ''nancolor'' is valid.',args{end-1})
    end
end

end % function checkArgs


%% function lbmap
function map = lbmap(n,scheme)
%LBMAP Returns specified Light-Bertlein colormap.
%
%   LBMAP(N,SCHEME) returns an Nx3 colormap. SCHEME can be one of the
%   following strings:
%
%       'Blue'       Single-hue progression to purlish-blue (default)
%       'BlueGray'   Diverging progression from blue to gray
%       'BrownBlue'  Orange-white-purple diverging scheme
%       'RedBlue'    Modified spectral scheme
%
%   If N is not specified, the size of the colormap is determined by the
%   current figure. If no figure exists, MATLAB creates one.
%
%Example 1: 7-color single-hue blue (default)
%   load penny
%   imagesc(P)
%   colormap(lbmap(7))
%   colorbar
%
%Example 2: 11-color modified spectrum
%   load penny
%   imagesc(P)
%   colormap(lbmap(11,'RedBlue'))
%   colorbar
%
%   See also HSV, GRAY, HOT, BONE, COPPER, PINK, FLAG, COLORMAP, RGBPLOT.

% Reference:
% A. Light & P.J. Bartlein, "The End of the Rainbow? Color Schemes for
% Improved Data Graphics," Eos,Vol. 85, No. 40, 5 October 2004.
% http://geography.uoregon.edu/datagraphics/EOS/Light&Bartlein_EOS2004.pdf

% Copyright 2007-2010 The MathWorks, Inc.

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
   
    
