function [] = save_data(data_dir,random_seed, rates_trial, ...
    times, tinds,params,stim_num)


file_name = sprintf('EIF_stim_num_%d_%d',stim_num,random_seed);
name_full = strcat(data_dir,file_name);
save(name_full,'times','tinds','rates_trial','params','-v7.3')

end

