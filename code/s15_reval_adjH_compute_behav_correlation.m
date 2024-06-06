% compute_behav_correlation.m

datapath = '/Users/nbertelsen/projects/cmi_eeg_rs_H/data';

rso_b1 = readtable(fullfile(datapath,'tidy','wide_rsopen_adjH_b1.csv'));
rso_b2 = readtable(fullfile(datapath,'tidy','wide_rsopen_adjH_b2.csv'));
rso_b3 = readtable(fullfile(datapath,'tidy','wide_rsopen_adjH_b3.csv'));
rso_b4 = readtable(fullfile(datapath,'tidy','wide_rsopen_adjH_b4.csv'));
rso_b5 = readtable(fullfile(datapath,'tidy','wide_rsopen_adjH_b5.csv'));

rso = [table2array(rso_b1(:,2:end)), ...
    table2array(rso_b2(:,2:end)), ...
    table2array(rso_b3(:,2:end)), ...
    table2array(rso_b4(:,2:end)), ...
    table2array(rso_b5(:,2:end))];

[coeff, rso_scores] = pca(rso);
rso_pc1 = [rso_b1.Var1, array2table(rso_scores(:,1))];
rso_pc1.Properties.RowNames = rso_pc1.Var1_1;
rso_pc1.Properties.VariableNames = {'subid','rso_PC1'};

rsc_b1 = readtable(fullfile(datapath,'tidy','rsc_b1.csv'));
rsc_b2 = readtable(fullfile(datapath,'tidy','rsc_b2.csv'));
rsc_b3 = readtable(fullfile(datapath,'tidy','rsc_b3.csv'));
rsc_b4 = readtable(fullfile(datapath,'tidy','rsc_b4.csv'));
rsc_b5 = readtable(fullfile(datapath,'tidy','rsc_b5.csv'));

rsc = [table2array(rsc_b1(:,2:end)), ...
    table2array(rsc_b2(:,2:end)), ... 
    table2array(rsc_b3(:,2:end)), ...
    table2array(rsc_b4(:,2:end)), ...
    table2array(rsc_b5(:,2:end))];

[coeff, rsc_scores] = pca(rsc);
rsc_pc1 = [rsc_b1.Var1, array2table(rsc_scores(:,1))];
rsc_pc1.Properties.RowNames = rsc_pc1.Var1_1;
rsc_pc1.Properties.VariableNames = {'subid','rsc_PC1'};


subs2use = rso_b1.Var1;

pheno_data = readtable(fullfile(datapath,'raw','rsopen_pheno.csv'));
pheno_data_subset = pheno_data(ismember(pheno_data.subid,subs2use), ...
    {'subid','subtype','age','srs_total_raw', ...
    'srs_socialawareness_T', 'srs_soccommint_T', ...
    'srs_socialcognition_T','srs_socialmotivation_T', ...
    'srs_socialcommunication_T', ...
    'srs_rrb_T',...
    'srs_total_T','scq_total'});
pheno_data_subset.Properties.RowNames = pheno_data_subset.subid;
pheno_data_subset = join(pheno_data_subset, rso_pc1);
pheno_data_subset = join(pheno_data_subset, rsc_pc1);

vars2use = {'rso_PC1','rsc_PC1','scq_total','srs_total_raw','srs_total_T','srs_soccommint_T','srs_rrb_T'};
[r_td,p_td] = corr(table2array(pheno_data_subset(ismember(pheno_data_subset.subtype,'td'), ...
    vars2use)));

[r_a1,p_a1] = corr(table2array(pheno_data_subset(ismember(pheno_data_subset.subtype,'autism1'), ...
    vars2use)));

[r_a2,p_a2] = corr(table2array(pheno_data_subset(ismember(pheno_data_subset.subtype,'autism2'), ...
    vars2use)));

[r_td(1,:);r_a1(1,:);r_a2(1,:)]
[r_td(2,:);r_a1(2,:);r_a2(2,:)]