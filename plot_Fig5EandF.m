%%
% Plot the firing rate and coherence curves for the VIP -> PV model
%
% Corresponds to Fig. 5E and F in:
%   J Veit, G Handy, DP Mossing, B Doiron, H Adesnik. 
%   Cortical VIP neurons locally control the gain but globally control 
%   the coherence of gamma band rhythms.
%
% Written by Gregory Handy, 08/24/2021
%%

clear; close all; clc;

%%
load('./Data_sets/EIF_LR_MF_data_vip_pv_surround_diff_0.mat');

%% Plot the firing rates (noting the oscillatory regime)
color_scheme = [38 30 101; 0 169 69; 255,172,16;...
    38 30 101; 0 169 69; 255,172,16]/255; 

figure(1); clf; hold on 
hh = plot(rates_trial_ave(1,1:7),'.-','markersize',16,'linewidth',1.5,'color',color_scheme(1,:));
gg = plot(8,[rates_min(1,8),rates_trial_ave(1,8)],'*','markersize',16,'linewidth',1.5,'color',color_scheme(1,:));

plot([7 8],[rates_trial_ave(1,7) rates_trial_ave(1,8)],'--','markersize',16,'linewidth',1.5,'color',color_scheme(1,:));
plot([7 8],[rates_trial_ave(1,7) rates_min(1,8)],'--','markersize',16,'linewidth',1.5,'color',color_scheme(1,:));

set(gca,'fontsize',16)
legend([hh gg(1)],{'stable','unstable'})
ylabel('Excitatory firing rates (Hz)')

xlabel('VIP Firing Rate')
xticks([1 8])
xticklabels({'Low', 'High'})

%% Plot the coherence
figure(2); clf; hold on;
% Find the peak of the normalized power curve
for ii = 1:7
    pos_indices = find(params.omega*1e3>14);
    [~, gamma_index(ii)] = max(real(squeeze(yy_freq(1,1,pos_indices,ii))));    
end

color_start = [0.949 0.8784 0.9451];
color_end = [0.6157 0.0549 0.5765];
color_interp = [];
for ii = 1:3
    color_interp = [color_interp; color_start(ii):...
        (color_end(ii)-color_start(ii))/7:color_end(ii)];
end
color_interp = color_interp';

subplot(1,2,1); hold on
for ii = 1:7    
    numerator = abs(squeeze(yy_freq(1,4,:,ii))).^2;
    denominator = real(squeeze(yy_freq(1,1,:,ii))).*real(squeeze(yy_freq(4,4,:,ii)));
    
    plot(params.omega*1e3,numerator./denominator,...
        '-','linewidth',1.5,'color',[0.6157 0.0549 0.5765 ii/8])
    xlim([0 100])
    ylim([0 1])
    set(gca,'fontsize',16)
    xlabel('Freq (Hz)')
    ylabel('coherence')
    
    coherence = numerator./denominator;
    
    save_max_co(ii) = coherence(pos_indices(gamma_index(ii)));
end

subplot(1,2,2); hold on
for ii = 1:7
    plot(ii,save_max_co(ii),'.','markersize',16,'color',color_interp(ii,:))
    
    if ii < 7
        plot([ii ii+1],[save_max_co(ii) save_max_co(ii+1)],'-','linewidth',1.5,'color',color_interp(ii,:))
    end
end

set(gca,'fontsize',16)
xlabel('Mod number')
ylabel('Gamma coherence')

ylim([0.7 1.0])
xlabel('VIP Firing Rate')
set(gca,'fontsize',16)
xticks([1 7])
xticklabels({'Low','High'})
xlim([1 7])