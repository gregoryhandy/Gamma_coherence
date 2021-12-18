function [] = save_data_corr_data(mod_num,random_seed,lags,corrAO_sim,...
    wfsp,sfsp)


file_name = sprintf('./Data_sets/Sim_%d/EIF_stim_num_%d_%d_cross_corr',...
    random_seed,mod_num,random_seed);
save(file_name,'lags','corrAO_sim','wfsp','sfsp')

end

