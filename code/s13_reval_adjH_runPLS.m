%% Run the PLS analysis
%ANALYSIS = 'rsclosed';

runPLS('rsopen')
runPLS('rsclosed')

function runPLS(ANALYSIS)
% runPLS
%
% ANALYSIS = 'rso' or 'rsc'

addpath /Users/nbertelsen/Documents/MATLAB/spm12;
addpath /Users/nbertelsen/Documents/MATLAB/pls/plsgui;
addpath /Users/nbertelsen/Documents/MATLAB/pls/plscmd;

rootpath = '/Users/nbertelsen/projects/cmi_eeg_rs_H';
codepath = fullfile(rootpath,'code');
datapath = fullfile(rootpath,'data','tidy');
resultpath = fullfile(rootpath,'results/reval/global');

RUN_BLOCKS = 0;

%% Prepare input arguments for pls_analysis.m
cd(codepath);

nblocks = 5;
nelectrodes = 93;

big_concat_data = cell(1,3);
big_concat_data{1} = []; 
big_concat_data{2} = []; 
big_concat_data{3} = []; 

tmp_data = readtable(fullfile(datapath, sprintf('wide_%s_adjH_b%d_td.csv',ANALYSIS,1)));
electrode_names = tmp_data.Properties.VariableNames(2:end);
electrode_labels = repmat(electrode_names', nblocks,1);
block_labels = [repmat({'block 1'}, length(electrode_names), 1); ...
    repmat({'block 2'}, length(electrode_names), 1); ...
    repmat({'block 3'}, length(electrode_names), 1); ...
    repmat({'block 4'}, length(electrode_names), 1); ...
    repmat({'block 5'}, length(electrode_names), 1)]; 


for iblock = 1:nblocks

    td_data = readtable(fullfile(datapath, sprintf('wide_%s_adjH_b%d_td.csv',ANALYSIS,iblock)));
    autism1_data = readtable(fullfile(datapath, sprintf('wide_%s_adjH_b%d_autism1.csv',ANALYSIS,iblock)));
    autism2_data = readtable(fullfile(datapath, sprintf('wide_%s_adjH_b%d_autism2.csv',ANALYSIS,iblock)));
    
    datamat_lst{iblock,1} = table2array(td_data(:,2:end));
    datamat_lst{iblock,2} = table2array(autism1_data(:,2:end));
    datamat_lst{iblock,3} = table2array(autism2_data(:,2:end));

    big_concat_data{1} = [big_concat_data{1} table2array(td_data(:,2:end))]; 
    big_concat_data{2} = [big_concat_data{2} table2array(autism1_data(:,2:end))]; 
    big_concat_data{3} = [big_concat_data{3} table2array(autism2_data(:,2:end))]; 

end % for iblock

% number of subjects per stacked data mat
num_subj_lst = [size(datamat_lst{1,1},1), size(datamat_lst{1,2},1), size(datamat_lst{1,3},1)];

pheno_fn = sprintf('tidy_%s_pheno_4pls.csv', ANALYSIS);
pheno_data = readtable(fullfile(datapath,pheno_fn));
pheno_mat = table2array(pheno_data(:,3:end));
pheno_names = pheno_data.Properties.VariableNames(3:end);

pheno_mat = zscore(pheno_mat);


num_cond = 1; %{ones(1,length(datamat_lst))};  % initialize num_cond

% specify option structure
% option.progress_hdl = [];%( user interface handle )
option.method = 3; %[1] | 2 | 3 | 4 | 5 | 6
option.num_perm = 10000; %( single non-negative integer )
option.is_struct = 0;%[0] | 1
option.num_split = 0; %( single non-negative integer )
option.num_boot = 10000; %( single non-negative integer )
option.clim = 95; %( [95] single number between 0 and 100 )
% option.bscan = ( subset of  1:num_cond )
% option.stacked_designdata = ( 2-D numerical matrix )
option.stacked_behavdata = pheno_mat;
% option.meancentering_type = [0] | 1 | 2 | 3
option.cormode = 0; %[0] | 2 | 4 | 6
option.boot_type = 'strat'; %['strat'] | 'nonstrat'


%% run pls_analysis.m
if RUN_BLOCKS

    for iblock = 1:nblocks

        % run analysis
        data2use = datamat_lst(iblock,:);
        result = pls_analysis(data2use, num_subj_lst, num_cond, option);

        % compute percentage of cross-block covariance
        result.crossblockCovPercent = result.s.^2/sum(result.s.^2);

        % compute the correct p-values
        result.perm_result.sprob_correct = (result.perm_result.sp+1)/(result.perm_result.num_perm+1);

        % find the significant LVs
        result.sigLVs = find(result.perm_result.sprob_correct<=0.05);
        result.sigLVs_pvals = result.perm_result.sprob_correct(result.sigLVs);
        
        disp(sprintf('%s block %d',ANALYSIS, iblock));
        disp(sprintf('Significant LVs: LV %d',result.sigLVs));
        disp(sprintf('Percentage of crossblock covariance explained = %f',result.crossblockCovPercent(result.sigLVs)));

        % compute brain BSR
        for i = 1:size(result.boot_result.compare_u,2)
            result.brain_bsr(:,i) = result.boot_result.compare_u(:,i)./result.boot_result.u_se(:,i);
        end

        % save result
        fname = fullfile(resultpath,sprintf('13_pls_%s_b%d.mat',ANALYSIS,iblock));
        save(fname,'result');

        % compute bootstrap CIs
        cd(codepath);
        LVnum = result.sigLVs;
        resultmat = fname;
        [tab2write] = s14_reval_adjH_compute_bootstrap_ci_pls(resultmat, LVnum, ANALYSIS, iblock, pheno_names);

    end % for iblock

end % if RUN_BLOCKS



result = pls_analysis(big_concat_data, num_subj_lst, num_cond, option);

% compute percentage of cross-block covariance
result.crossblockCovPercent = result.s.^2/sum(result.s.^2);

% compute the correct p-values
result.perm_result.sprob_correct = (result.perm_result.sp+1)/(result.perm_result.num_perm+1);

% find the significant LVs
result.sigLVs = find(result.perm_result.sprob_correct<=0.05);
result.sigLVs_pvals = result.perm_result.sprob_correct(result.sigLVs);

disp(sprintf('%s block %d',ANALYSIS, iblock));
disp(sprintf('Significant LVs: LV %d',result.sigLVs));
disp(sprintf('Percentage of crossblock covariance explained = %f',result.crossblockCovPercent(result.sigLVs)));

% compute brain BSR
for i = 1:size(result.boot_result.compare_u,2)
    result.brain_bsr(:,i) = result.boot_result.compare_u(:,i)./result.boot_result.u_se(:,i);
end

tab2write = cell2table([electrode_labels, ...
    block_labels, ...
    num2cell(result.brain_bsr(:,result.sigLVs))], ...
    'VariableNames',{'electrode','block','BSR'});
bsrname2save = fullfile(resultpath, sprintf('13_pls_%s_ALL_BSR_LV%d.csv',ANALYSIS,result.sigLVs));
writetable(tab2write,bsrname2save);

tab2write = cell2table([electrode_labels, ...
    block_labels, ...
    num2cell(result.brain_bsr(:,result.sigLVs).*-1)], ...
    'VariableNames',{'electrode','block','BSR'});
bsrname2save = fullfile(resultpath, sprintf('13_pls_%s_ALL_BSRrev_LV%d.csv',ANALYSIS,result.sigLVs));
writetable(tab2write,bsrname2save);

% save result
fname = fullfile(resultpath,sprintf('13_pls_%s_ALL.mat',ANALYSIS));
save(fname,'result');

% compute bootstrap CIs
cd(codepath);
LVnum = result.sigLVs;
resultmat = fname;
[tab2write] = s14_reval_adjH_compute_bootstrap_ci_pls(resultmat, LVnum, ANALYSIS, 'ALL', pheno_names);


end % function runPLS
