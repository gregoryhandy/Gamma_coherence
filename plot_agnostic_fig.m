clear; close all;clc;


cross_theory = load('./Data_sets/EIF_LR_MF_data_vip_all_surround_diff_90.mat');
iso_theory = load('./Data_sets/EIF_LR_MF_data_vip_all_surround_diff_0.mat');

%%
color_start = [0.949 0.8784 0.9451];
color_end = [0.6157 0.0549 0.5765];


color_interp = [];
for ii = 1:3
color_interp = [color_interp; color_start(ii):(color_end(ii)-color_start(ii))/7:color_end(ii)];
end
color_interp = color_interp';


indices_to_plot = [1:8];
figure(33); clf; 
subplot(1,2,1); hold on
for kk = 1:length(indices_to_plot)
    ii = indices_to_plot(kk);
    
    numerator = abs(squeeze(cross_theory.yy_freq(1,4,:,ii))).^2;
    denominator = real(squeeze(cross_theory.yy_freq(1,1,:,ii))).*real(squeeze(cross_theory.yy_freq(4,4,:,ii)));
    
    if kk == 2 || kk == 8
    h_blah(kk) = plot(cross_theory.params.omega*1e3,numerator./denominator,...
        '-','linewidth',1.5,'color',[0.6157 0.0549 0.5765 ii/8]);
    end
    xlim([0 100])
    ylim([0 1])
    set(gca,'fontsize',16)
    xlabel('Freq (Hz)')
    ylabel('coherence')
    
    if strcmp('Cross','Iso')
        gamma_range = [find(cross_theory.params.omega*1e3>0,1) find(cross_theory.params.omega*1e3>100,1)];
    else
        gamma_range = [find(cross_theory.params.omega*1e3>25,1) find(cross_theory.params.omega*1e3>60,1)];
    end
    
    cohere_range = (numerator(gamma_range(1):gamma_range(2))./...
        denominator(gamma_range(1):gamma_range(2)));
    freq_range = cross_theory.params.omega(gamma_range(1):gamma_range(2))*1e3;
    
    [save_max_co_cross(kk),gamma_index(kk)] = max(cohere_range);
   
end


% subplot(2,2,4); hold on
% for ii = 1:length(indices_to_plot)
%     plot(ii,save_max_co_cross(ii),'.','markersize',16,'color',color_interp(ii,:))
%     
%     if ii < length(indices_to_plot)
%         plot([ii ii+1],[save_max_co_cross(ii) save_max_co_cross(ii+1)],'-','linewidth',1.5,'color',color_interp(ii,:))
%     end
% end
% 
% set(gca,'fontsize',16)
% xlabel('Mod number')
% ylabel('Gamma coherence')
% 
% xlabel('VIP Firing Rate')
% set(gca,'fontsize',16)
% xticks([1 indices_to_plot(end)])
% xticklabels({'Low','High'})
% xlim([1 indices_to_plot(end)])
% 
% sgtitle(sprintf('%s Surround','Cross'),'FontSize',16,'FontWeight','bold')



indices_to_plot = [1:8];
for kk = 1:length(indices_to_plot)
    ii = indices_to_plot(kk);
    
    numerator = abs(squeeze(iso_theory.yy_freq(1,4,:,ii))).^2;
    denominator = real(squeeze(iso_theory.yy_freq(1,1,:,ii))).*real(squeeze(iso_theory.yy_freq(4,4,:,ii)));
    
    if kk == 2 || kk == 8
        plot(iso_theory.params.omega*1e3,numerator./denominator,...
        '--','linewidth',1.5,'color',[0.6157 0.0549 0.5765 ii/8])
    end
    xlim([0 100])
    ylim([0 1])
    set(gca,'fontsize',16)
    xlabel('Freq (Hz)')
    ylabel('coherence')
    
    if strcmp('Iso','Iso')
        gamma_range = [find(iso_theory.params.omega*1e3>0,1) find(iso_theory.params.omega*1e3>100,1)];
    else
        gamma_range = [find(iso_theory.params.omega*1e3>25,1) find(iso_theory.params.omega*1e3>60,1)];
    end
    
    cohere_range = (numerator(gamma_range(1):gamma_range(2))./...
        denominator(gamma_range(1):gamma_range(2)));
    freq_range = cross_theory.params.omega(gamma_range(1):gamma_range(2))*1e3;
    
    [save_max_co_iso(kk),gamma_index(kk)] = max(cohere_range);
    
end

legend([h_blah(2),h_blah(8)],{'VIP Low','VIP High'}) 

subplot(1,2,2); hold on
for ii = 1:length(indices_to_plot)
    plot(ii,save_max_co_iso(ii),'.','markersize',16,'color',color_interp(ii,:))
    
    if ii < length(indices_to_plot)
         h_iso(ii) = plot([ii ii+1],[save_max_co_iso(ii) save_max_co_iso(ii+1)],'--','linewidth',1.5,'color',color_interp(ii,:));
    end
end

for ii = 1:length(indices_to_plot)
    plot(ii,save_max_co_cross(ii),'.','markersize',16,'color',color_interp(ii,:))
    
    if ii < length(indices_to_plot)
        h_cross(ii) = plot([ii ii+1],[save_max_co_cross(ii) save_max_co_cross(ii+1)],'-','linewidth',1.5,'color',color_interp(ii,:));
    end
end

ylabel('Gamma coherence')
xlabel('VIP Firing Rate')
set(gca,'fontsize',16)
xticks([1 indices_to_plot(end)])
xticklabels({'Low','High'})
xlim([1 indices_to_plot(end)])
legend([h_iso(end) h_cross(end)], {'Iso','Cross'})
