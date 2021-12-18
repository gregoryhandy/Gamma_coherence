clear; close all; clc;

surround_diff = 90; % 0 (iso), 90 (cross) or 180 (2-pop)
if surround_diff == 0
    theory_results = load('./Data_sets/EIF_LR_MF_data_vip_surround_diff_0.mat');
elseif surround_diff == 90
    theory_results = load('./Data_sets/EIF_LR_MF_data_vip_surround_diff_90.mat');
else
    theory_results = load('./Data_sets/EIF_LR_MF_data_2pop_vip.mat');
end
params = theory_results.params;

%%

color_scheme = [0.8500, 0.3250, 0.0980; 0, 0.4470, 0.7410; 0.4660, 0.6740, 0.1880;...
    0.8500, 0.3250, 0.0980; 0, 0.4470, 0.7410; 0.4660, 0.6740, 0.1880]; 

figure(1); clf;
for ii = 1:theory_results.params.Npop
    h(ii) = plot(theory_results.rates_trial_ave(ii,:),'.-','markersize',16,'linewidth',1.5,'color',color_scheme(ii,:));
    hold on
    
    if exist('condensed_sim1','var')
        plot(sim_rates(ii,:),'*','markersize',10,'linewidth',1.5,'color',color_scheme(ii,:))
    end
    
    if exist('LR_theory','var')
        plot(LR_theory.rates_trial_ave(ii,:),'o','markersize',10,'linewidth',1.5,'color',color_scheme(ii,:))
    end
end
set(gca,'fontsize',16)
legend(h(1:3),{'E','PV','SOM'})
ylabel('Firing rate (Hz)')


%%
figure(22); clf; hold on;
indices_to_plot = [1:8];
for kk = 1:length(indices_to_plot)
    ii = indices_to_plot(kk);
    norm_factor = 1+0*real(squeeze(theory_results.yy_freq(1,1,params.ind0,ii)));
    plot(params.omega*1e3,(real(squeeze(theory_results.yy_freq(1,1,:,ii)))/norm_factor),...
        '-','linewidth',1.5,'color',[0 0.4470 0.7410 ii/8])
    xlim([0 100])
    
    pos_indices = find(params.omega*1e3>0);
    [val, index] = max(real(squeeze(theory_results.yy_freq(1,1,pos_indices,ii))));
    max_freq(ii) = params.omega(pos_indices(index))*1e3;
end
set(gca,'fontsize',16)
xlabel('Freq (Hz)')
ylabel('c_{e1e1}(f)')


%%
clearvars h;
figure(2); clf; hold on;
subplot(2,2,1); hold on

% indices_to_plot = [1:8];
indices_to_plot = [1:7];
for kk = 1:length(indices_to_plot)
    ii = indices_to_plot(kk);
    norm_factor = real(squeeze(theory_results.yy_freq(1,1,params.ind0,ii)));
    h(kk) = plot(params.omega*1e3,(real(squeeze(theory_results.yy_freq(1,1,:,ii)))/norm_factor),...
        '-','linewidth',1.5,'color',[0.6157 0.0549 0.5765 ii/length(indices_to_plot)]);
    xlim([0 100])
    
    pos_indices = find(params.omega*1e3>14);
    [gamma_power(ii), gamma_index(ii)] = max(real(squeeze(theory_results.yy_freq(1,1,pos_indices,ii)))/norm_factor);
    max_freq(ii) = params.omega(pos_indices(index))*1e3;
    
end
set(gca,'fontsize',16)
xlabel('Freq (Hz)')
ylabel('Power (normalized)')
legend([h(1),h(end)],{'VIP Low','VIP High'}) 

subplot(2,2,2); hold on
color_start = [0.949 0.8784 0.9451];
color_end = [0.6157 0.0549 0.5765];
color_interp = [];
for ii = 1:3
color_interp = [color_interp; color_start(ii):(color_end(ii)-color_start(ii))/(length(indices_to_plot)-1):color_end(ii)];
end
color_interp = color_interp';
for ii = 1:length(indices_to_plot)
    plot(ii,gamma_power(ii),'.','markersize',16,'color',color_interp(ii,:))
    
    if ii < length(indices_to_plot)
        plot([ii ii+1],[gamma_power(ii) gamma_power(ii+1)],'-','linewidth',1.5,'color',color_interp(ii,:))
    end
end
ylabel('Normalized gamma power')
xlabel('VIP Firing Rate')
set(gca,'fontsize',16)
xticks([indices_to_plot(1) indices_to_plot(end)])
xticklabels({'Low','High'})
xlim([indices_to_plot(1) indices_to_plot(end)])


subplot(2,2,3); hold on
for kk = 1:length(indices_to_plot)
    ii = indices_to_plot(kk);
    
    numerator = abs(squeeze(theory_results.yy_freq(1,4,:,ii))).^2;
    denominator = real(squeeze(theory_results.yy_freq(1,1,:,ii))).*real(squeeze(theory_results.yy_freq(4,4,:,ii)));
    
    plot(params.omega*1e3,numerator./denominator,...
        '-','linewidth',1.5,'color',[0.6157 0.0549 0.5765 ii/length(indices_to_plot)])
    xlim([0 100])
    ylim([0 1])
    set(gca,'fontsize',16)
    xlabel('Freq (Hz)')
    ylabel('coherence')
    
%     gamma_range = [find(params.omega*1e3>14,1) find(params.omega*1e3>100,1)];
%     cohere_range = (numerator(gamma_range(1):gamma_range(2))./...
%         denominator(gamma_range(1):gamma_range(2)));
%     freq_range = params.omega(gamma_range(1):gamma_range(2))*1e3;

    coherence = numerator./denominator;
    [save_max_co(kk)] = coherence(pos_indices(gamma_index(ii)));

%     [save_max_co(kk),gamma_index(kk)] = max(cohere_range);
end

subplot(2,2,4); hold on
for ii = 1:length(indices_to_plot)
    plot(ii,save_max_co(ii),'.','markersize',16,'color',color_interp(ii,:))
    
    if ii < length(indices_to_plot)
        plot([ii ii+1],[save_max_co(ii) save_max_co(ii+1)],'-','linewidth',1.5,'color',color_interp(ii,:))
    end
end

set(gca,'fontsize',16)
xlabel('Mod number')
ylabel('Gamma coherence')

if surround_diff==0
    ylim([0.35 0.65])
elseif surround_diff == 90
    ylim([0.2 0.5])
end

xlabel('VIP Firing Rate')
set(gca,'fontsize',16)
xticks([1 indices_to_plot(end)])
xticklabels({'Low','High'})
xlim([1 indices_to_plot(end)])

