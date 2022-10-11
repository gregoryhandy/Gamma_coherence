%% 
% Compares and plots the theory vs. simulation (using simulation data
% already acquired)
%
% Corresponds to Fig. 4C and D in:
%   J Veit, G Handy, DP Mossing, B Doiron, H Adesnik. 
%   Cortical VIP neurons locally control the gain but globally control 
%   the coherence of gamma band rhythms.
% Also produces additional matches for other VIP modulation values
%
% Written by Gregory Handy, 08/24/2021
%%
clear; clc; close all;

%% Load the simulation data and the mean field data
random_seed = 1635400248; % iso surround simulations
surround_diff = 0;

LF_mf_file_name = sprintf('./Data_sets/EIF_LR_MF_data_vip_surround_diff_%d.mat',surround_diff);
LR_mean_field_data = load(LF_mf_file_name);

sim_data_condensed = load(sprintf('./Data_sets/Sim_%d/EIF_condensed_%d.mat',random_seed,random_seed));

%%
color_scheme = [38 30 101; 0 169 69; 255,172,16;...
    38 30 101; 0 169 69; 255,172,16]/255; 

figure(1); clf; 
subplot(1,2,1); hold on;
for ii = 1:3
    h(ii) = plot(LR_mean_field_data.rates_trial_ave(ii,:),'linewidth',3,'color',color_scheme(ii,:));
    plot(LR_mean_field_data.rates_trial_ave(ii+3,:),'linewidth',3,'color',color_scheme(ii,:))
    g(ii) = plot(sim_data_condensed.rates_ave(ii,:),'.','markersize',20,'color',color_scheme(ii,:));
    plot(sim_data_condensed.rates_ave(ii+3,:),'.','markersize',20,'color',color_scheme(ii,:))
end
set(gca,'fontsize',30)
ylabel('Firing rates')
xlabel('VIP Firing Rate')
xticks([1 8])
xlim([1 8])
xticklabels({'Low', 'High'})
legend([h(1) g(1)],'Mean Field Theory','Spiking Simulation')



for ii = 5:5
    mod_num1 = ii;
    file_name = sprintf('./Data_sets/Sim_%d/EIF_stim_num_%d_%d_power_spec.mat',random_seed,mod_num1,random_seed);
    sim_data = load(file_name);
    
    subplot(1,2,2); hold on;
    plot(LR_mean_field_data.params.omega*1e3,real(squeeze(LR_mean_field_data.yy_freq(1,1,:,mod_num1))),...
        'k-','linewidth',1.5);
    plot(sim_data.wfsp1/(2*pi)*1e3,sim_data.sfsp1,'k.','markersize',10)    
    xlim([0 100])
    if ii >12
        xlabel('Freq (Hz)')
    end
    if mod(ii,4)==1 
        ylabel('c_{ee}(f)')
    end
    ylabel('c_{ee}(f) (spikes^2 per Hz)')
    xlabel('Freq (Hz)')
    set(gca,'fontsize',30)
    if ii == 1
        legend([h(1) g(1)],{'Theory','Sim'})
    end
end

%%

figure(3); clf;
for ii = 1:8
    
    mod_num1 = ii;
    file_name = sprintf('./Data_sets/Sim_%d/EIF_stim_num_%d_%d_cross_spec.mat',random_seed,mod_num1,random_seed);
    sim1_cross_spec = load(file_name);
    
    file_name = sprintf('./Data_sets/Sim_%d/EIF_stim_num_%d_%d_power_spec.mat',random_seed,mod_num1,random_seed);
    sim1_power_spec = load(file_name);
    
    theory_coherence = abs(squeeze(LR_mean_field_data.yy_freq(1,4,:,mod_num1))).^2./...
        (real(squeeze(LR_mean_field_data.yy_freq(1,1,:,mod_num1))).*...
        real(squeeze(LR_mean_field_data.yy_freq(4,4,:,mod_num1))));
        
    
    sim_coherence = abs(sim1_cross_spec.sfsp).^2./(sim1_power_spec.sfsp1.*sim1_power_spec.sfsp2);
    
    subplot(2,4,ii); hold on;
    plot(LR_mean_field_data.params.omega*1e3,theory_coherence,...
        '-','linewidth',1.5)
    plot(sim1_cross_spec.wfsp/(2*pi)*1e3,sim_coherence,'.','linewidth',1.5)
    
    if ii == 1
        gamma_index_theory = find(LR_mean_field_data.params.omega*1e3>12,1);
        gamma_index_sim = find(sim1_cross_spec.wfsp/(2*pi)*1e3>12,1);
    end
    
    gamma_cohere_theory(ii) = theory_coherence(gamma_index_theory);
    gamma_cohere_sim(ii) = sim_coherence(gamma_index_sim);
    
    xlim([0 100])
    if ii >12
        xlabel('Freq (Hz)')
    end
    if mod(ii,4)==1 
        ylabel('coherence')
    end
    set(gca,'fontsize',16)
    if ii == 1
        legend('Theory','Sim')
    end
    
    title(sprintf('Mod num %d',ii))
end


