%% cd Theta
rho_v = ['0.0' string(0.1:0.1:0.9) '1.0'];
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

writetable(stats_all, fullfile('outputs', 'stats_all.csv'), 'WriteRowNames', false);
