%%
% Plot the normalized power and coherence curves for different surround
% orientations using the theory
%
% Corresponds to Fig. 5A-D in:
%   J Veit, G Handy, DP Mossing, B Doiron, H Adesnik. 
%   Cortical VIP neurons locally control the gain but globally control 
%   the coherence of gamma band rhythms.
%
% Written by Gregory Handy, 08/24/2021
%%
clear; close all; clc;

% Loop through the two conditions and plot the results
surround_diffs = [0 90]; % 0 (iso), 90 (cross)
for outerLoop = 1:2
    
    surround_diff = surround_diffs(outerLoop);
    if surround_diff == 0
        theory_results = load('./Data_sets/EIF_LR_MF_data_vip_surround_diff_0.mat');
    elseif surround_diff == 90
        theory_results = load('./Data_sets/EIF_LR_MF_data_vip_surround_diff_90.mat');
    else
        error('Something went wrong loading the data.')
    end
    params = theory_results.params;
    
    %%
    
    color_scheme = [0.8500, 0.3250, 0.0980; 0, 0.4470, 0.7410; 0.4660, 0.6740, 0.1880;...
        0.8500, 0.3250, 0.0980; 0, 0.4470, 0.7410; 0.4660, 0.6740, 0.1880];
    
    %%
    figure(outerLoop); clf;
    
    %% Plot the normalized power curves
    subplot(2,2,1); hold on
    indices_to_plot = [1:7];
    for kk = 1:length(indices_to_plot)
        ii = indices_to_plot(kk);
        norm_factor = real(squeeze(theory_results.yy_freq(1,1,params.ind0,ii)));
        h(kk) = plot(params.omega*1e3,(real(squeeze(theory_results.yy_freq(1,1,:,ii)))/norm_factor),...
            '-','linewidth',1.5,'color',[0.6157 0.0549 0.5765 ii/length(indices_to_plot)]);
        xlim([0 100])
        
        pos_indices = find(params.omega*1e3>14);
        [gamma_power(ii), gamma_index(ii)] = max(real(squeeze(theory_results.yy_freq(1,1,pos_indices,ii)))/norm_factor);
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
        color_interp = [color_interp; color_start(ii):...
            (color_end(ii)-color_start(ii))/(length(indices_to_plot)-1):...
            color_end(ii)];
    end
    color_interp = color_interp';
    for ii = 1:length(indices_to_plot)
        plot(ii,gamma_power(ii),'.','markersize',16,'color',color_interp(ii,:))
        
        if ii < length(indices_to_plot)
            plot([ii ii+1],[gamma_power(ii) gamma_power(ii+1)],'-',...
                'linewidth',1.5,'color',color_interp(ii,:))
        end
    end
    ylabel('Normalized gamma power')
    xlabel('VIP Firing Rate')
    set(gca,'fontsize',16)
    xticks([indices_to_plot(1) indices_to_plot(end)])
    xticklabels({'Low','High'})
    xlim([indices_to_plot(1) indices_to_plot(end)])
    
    %% Plot the coherence curves
    subplot(2,2,3); hold on
    for kk = 1:length(indices_to_plot)
        ii = indices_to_plot(kk);
        
        numerator = abs(squeeze(theory_results.yy_freq(1,4,:,ii))).^2;
        denominator = real(squeeze(theory_results.yy_freq(1,1,:,ii))).*...
            real(squeeze(theory_results.yy_freq(4,4,:,ii)));
        
        plot(params.omega*1e3,numerator./denominator,...
            '-','linewidth',1.5,'color',[0.6157 0.0549 0.5765 ii/length(indices_to_plot)])
        xlim([0 100])
        ylim([0 1])
        set(gca,'fontsize',16)
        xlabel('Freq (Hz)')
        ylabel('coherence')
        
        coherence = numerator./denominator;
        [save_max_co(kk)] = coherence(pos_indices(gamma_index(ii)));
    end
    
    subplot(2,2,4); hold on
    for ii = 1:length(indices_to_plot)
        plot(ii,save_max_co(ii),'.','markersize',16,'color',color_interp(ii,:))
        
        if ii < length(indices_to_plot)
            plot([ii ii+1],[save_max_co(ii) save_max_co(ii+1)],'-','linewidth',1.5,'color',color_interp(ii,:))
        end
    end
    
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
    ylabel('Gamma coherence')
end

