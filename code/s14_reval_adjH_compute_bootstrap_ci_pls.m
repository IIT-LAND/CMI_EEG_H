function [tab2write] = s14_reval_adjH_compute_bootstrap_ci_pls(resultmat, LVnum, ANALYSIS, iblock, pheno_names)
%
%   resultmat = full path and filename to PLS resultmat
%   LVnum = number of LV you want
%

ORDERBY = 'TD';
mod_names_orig = pheno_names;

mod_names = repmat(mod_names_orig',3,1);

[fpath, fname, fext] = fileparts(resultmat);
load(resultmat);

td_idx = 1:length(mod_names_orig);
autism1_idx = (length(mod_names_orig)+1):(length(mod_names_orig)*2);
autism2_idx = ((length(mod_names_orig)*2)+1):(length(mod_names_orig)*3);

td_bootres = squeeze(result.boot_result.distrib(td_idx,LVnum,:));
autism1_bootres = squeeze(result.boot_result.distrib(autism1_idx,LVnum,:));
autism2_bootres = squeeze(result.boot_result.distrib(autism2_idx,LVnum,:));

ci_bounds = [2.5,97.5];
autism1_ci = prctile(autism1_bootres',ci_bounds)';
autism2_ci = prctile(autism2_bootres',ci_bounds)';
td_ci = prctile(td_bootres',ci_bounds)';

autism1_corr = result.boot_result.orig_corr(autism1_idx,LVnum);
autism2_corr = result.boot_result.orig_corr(autism2_idx,LVnum);
td_corr = result.boot_result.orig_corr(td_idx,LVnum);

if strcmp(ORDERBY,'TD')
    [idx, plot_order] = sort(td_corr,'ascend');
    plot_order(plot_order) = 1:length(mod_names_orig);
elseif strcmp(ORDERBY,'Autism1')
    [idx, plot_order] = sort(autism1_corr,'ascend');
    plot_order(plot_order) = 1:length(mod_names_orig);
elseif strcmp(ORDERBY,'Autism2')
    [idx, plot_order] = sort(autism2_corr,'ascend');
    plot_order(plot_order) = 1:length(mod_names_orig);
end

autism1 = [autism1_corr autism1_ci];
autism2 = [autism2_corr autism2_ci];
td = [td_corr td_ci];

all_data = [td; autism1; autism2];
all_data(:,4) = (sign(all_data(:,1)) == sign(all_data(:,2))) & (sign(all_data(:,1)) == sign(all_data(:,3)));
all_data(:,5) = repmat(plot_order,3,1);


group_labels = [repmat({'TD'},length(mod_names_orig),1); ...
    repmat({'Autism 1'},length(mod_names_orig),1); ...
    repmat({'Autism 2'},length(mod_names_orig),1)];

tab2write = cell2table([group_labels, mod_names num2cell(all_data)], ...
    'VariableNames',{'Grp','VarName','corr','lo_lim','up_lim','nonzero','plot_order'});

if strcmp(iblock,'ALL')
    fname2save = fullfile(fpath,sprintf('14_%s_%s_bootlim_data4plotting_LV%d_ci%d.csv', ...
        ANALYSIS, iblock, LVnum, ci_bounds(2)-ci_bounds(1)));
else
    fname2save = fullfile(fpath,sprintf('14_%s_block%d_bootlim_data4plotting_LV%d_ci%d.csv', ...
        ANALYSIS, iblock, LVnum, ci_bounds(2)-ci_bounds(1)));
end % if strcmp(iblock,'ALL')
writetable(tab2write,fname2save,'FileType','text','delimiter',',');

all_data_rev = [td; autism1; autism2].*-1;
all_data_rev(:,4) = all_data(:,4);
all_data_rev(:,5) = all_data(:,5);

if strcmp(iblock,'ALL')
    
    tab2write = cell2table([group_labels, mod_names num2cell(all_data_rev)], ...
        'VariableNames',{'Grp','VarName','corr','lo_lim','up_lim','nonzero','plot_order'});
    
    fname2save = fullfile(fpath,sprintf('14_%s_%srev_bootlim_data4plotting_LV%d_ci%d.csv', ...
        ANALYSIS, iblock, LVnum, ci_bounds(2)-ci_bounds(1)));
    
    writetable(tab2write,fname2save,'FileType','text','delimiter',',');

end % if strcmp(iblock,'ALL')
