%% cd Theta

rho_v = ['0.0' string(0.1:0.1:0.9) '1.0'];
theta = sqrt(2)/10;
stats_all = [];
for rho = rho_v
    excrs = 'rho_'+rho+'_newtau';
    load('outputs/MC_'+excrs+'.mat','stats_tbl');
    stats_all = [ stats_all ; stats_tbl];
end
corr=str2double(rho_v).*exp(-(0.1/theta));
stats_all= addvars(stats_all, corr', 'Before', "Rej F",'NewVariableNames', 'Corr');
stats_all= addvars(stats_all, str2double(rho_v)', 'Before', "Corr",'NewVariableNames', '$\rho$');
stats_all

%% Save all stats

save(['outputs/collect_all_'+excrs+'.mat'])
writetable(stats_all, fullfile('outputs', 'stats_all.csv'), 'WriteRowNames', false);

%% Collecting PCs %%

rho_v = ['0.0' string(0.1:0.1:0.9) '1.0'];
pcs_all = [];
bws_all = [];
for rho = rho_v
    excrs = 'rho_'+rho+'_newtau';
    load('outputs/MC_'+excrs+'.mat','PC_n','MC_array','l_cutoffs');
    avg_PC = sum(MC_array,1,'omitnan')*PC_n'/100
    avg_bw = mean(MC_array,2,'omitnan')*l_cutoffs'/100
    pcs_all = [ pcs_all ; avg_PC];
    bws_all = [ bws_all ; avg_bw];
end

stats_all = addvars(stats_all, pcs_all, 'NewVariableNames', 'PCs');
stats_all
%% Save all stats

save(['outputs/collect_all_'+excrs+'.mat'])
writetable(stats_all, fullfile('outputs', 'stats_all_pcs.csv'), 'WriteRowNames', false);

%% Collecting CVs %%

rho_v = ['0.0' string(0.1:0.1:0.9) '1.0'];
cvs_all = [];
for rho = rho_v
    excrs = 'rho_'+rho+'_newtau';
    load('outputs/MC_'+excrs+'.mat','cv_array');
    avg_cv = mean(cv_array,'all','omitnan')
    cvs_all = [ cvs_all ; avg_cv];
end
load(['outputs/collect_all_'+excrs+'.mat'])
stats_all = addvars(stats_all, cvs_all, 'NewVariableNames', 'CVs');
stats_all

%% Save all stats

save(['outputs/collect_all_'+excrs+'.mat'])
writetable(stats_all, fullfile('outputs', 'stats_all_cvs.csv'), 'WriteRowNames', false);

%% Collecting BWs %%

rho_v = ['0.0' string(0.1:0.1:0.9) '1.0'];
% pcs_all = [];
bws_all = [];
for rho = rho_v
    excrs = 'rho_'+rho+'_newtau';
    load('outputs/MC_'+excrs+'.mat','PC_n','MC_array','l_cutoffs');
    % avg_PC = sum(MC_array,1,'omitnan')*PC_n'/100
    avg_bw = (sum(MC_array,2,'omitnan')'/100)*l_cutoffs'
    % pcs_all = [ pcs_all ; avg_PC];
    bws_all = [ bws_all ; avg_bw];
end
load(['outputs/collect_all_'+excrs+'.mat'],'stats_all')
stats_all = addvars(stats_all, bws_all, 'NewVariableNames', 'BWs');
stats_all
%% Save all stats

save(['outputs/collect_all_'+excrs+'.mat'])
writetable(stats_all, fullfile('outputs', 'stats_all_bws.csv'), 'WriteRowNames', false);