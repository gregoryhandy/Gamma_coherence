%%
% Creates the normalized power spectrum figures (using simulation data
% already acquired)
% Corresponds to Fig. 4E-J in:
%   J Veit, G Handy, DP Mossing, B Doiron, H Adesnik. 
%   Cortical VIP neurons locally control the gain but globally control 
%   the coherence of gamma band rhythms.
%
% Written by Gregory Handy, 08/24/2021
%%
clear; close all; clc;

%% Load the mean field data for all conditions
theory_size_control = load('./Data_sets/Exp_compare/EIF_LR_MF_data_size_light_off.mat');
theory_size_light = load('./Data_sets/Exp_compare/EIF_LR_MF_data_size_light_on.mat');

theory_contrast_control = load('./Data_sets/Exp_compare/EIF_LR_MF_data_contrast_light_off.mat');
theory_contrast_light = load('./Data_sets/Exp_compare/EIF_LR_MF_data_contrast_light_on.mat');

theory_surround_control = load('./Data_sets/Exp_compare/EIF_LR_MF_data_surround_light_off.mat');
theory_surround_light = load('./Data_sets/Exp_compare/EIF_LR_MF_data_surround_light_on.mat');

%%
clearvars gamma_power_control gamma_power_light hh gg

figure(11);
subplot(2,3,1); hold on;
indices_to_plot = [1:4];
for kk = 1:length(indices_to_plot)
    ii = indices_to_plot(kk);
    norm_factor = real(squeeze(theory_size_control.yy_freq(1,1,theory_size_control.params.ind0,ii)));
    
    if ii == 2 || ii == 4
        hh(ii) = plot(theory_size_control.params.omega*1e3,(real(squeeze(theory_size_control.yy_freq(1,1,:,ii)))/norm_factor),...
            '-','linewidth',1.5,'color',[0 0 0 ii/4]);
        xlim([0 100])
    end
    
    pos_indices = find(theory_size_control.params.omega*1e3>0);
    [gamma_power_control(ii), index] = max(real(squeeze(theory_size_control.yy_freq(1,1,pos_indices,ii)))/norm_factor);
    max_freq(ii) = theory_size_control.params.omega(pos_indices(index))*1e3;
end

indices_to_plot = [1:4];
for kk = 1:length(indices_to_plot)
    ii = indices_to_plot(kk);
    norm_factor = real(squeeze(theory_size_light.yy_freq(1,1,theory_size_light.params.ind0,ii)));
    
    if ii == 2 || ii == 4
        gg(kk) = plot(theory_size_light.params.omega*1e3,(real(squeeze(theory_size_light.yy_freq(1,1,:,ii)))/norm_factor),...
            '-','linewidth',1.5,'color',[1 0.0353 0 ii/4]);
        xlim([0 100])
    end
    
    pos_indices = find(theory_size_light.params.omega*1e3>0);
    [gamma_power_light(ii), index] = max(real(squeeze(theory_size_light.yy_freq(1,1,pos_indices,ii)))/norm_factor);
    max_freq(ii) = theory_size_light.params.omega(pos_indices(index))*1e3;
end
set(gca,'fontsize',16)
xlabel('Freq (Hz)')
ylabel('c_{e1e1}(f) (normalized)')
legend([hh(2) gg(2) hh(4) gg(4)],{'small','small (light)','large','large (light)'})

subplot(2,3,4); hold on
plot(gamma_power_control,'.-','markersize',16,'color','k',...
        'linewidth',1.5)
plot(gamma_power_light,'.-','markersize',16,'color',[1 0.0353 0],...
        'linewidth',1.5)
ylabel('Normalized gamma power')
xlabel('Stimulus Size')
set(gca,'fontsize',16)
xticks([1 indices_to_plot(2) 3 4])
xticklabels({'spon.', 'small','med','large'})
xlim([indices_to_plot(1) 4])
ylim_temp = ylim;
ylim([1 ylim_temp(2)])

figure(32); clf; 
subplot(1,3,1); hold on;
plot(gamma_power_control,gamma_power_light,'k.-','linewidth',2,'markersize',30)
x_temp = min(gamma_power_control):0.001:max(gamma_power_control);
plot(x_temp,x_temp,'k')
set(gca,'fontsize',16)
xlabel('Gamma Control')
ylabel('Gamma Light')
xlim([min(gamma_power_control) max(gamma_power_control)])
xticks([])
yticks([])

%%
clearvars gamma_power_control gamma_power_light hh gg
figure(11);
subplot(2,3,2); hold on;
indices_to_plot = [1:4];
for kk = 1:length(indices_to_plot)
    ii = indices_to_plot(kk);
    norm_factor = real(squeeze(theory_contrast_control.yy_freq(1,1,theory_contrast_control.params.ind0,ii)));
    
    if kk == 2 || kk == 4
        hh(kk) = plot(theory_contrast_control.params.omega*1e3,(real(squeeze(theory_contrast_control.yy_freq(1,1,:,ii)))/norm_factor),...
            '-','linewidth',1.5,'color',[0 0 0 ii/length(indices_to_plot)]);
        xlim([0 100])
    end
    
    pos_indices = find(theory_contrast_control.params.omega*1e3>0);
    [gamma_power_control(ii), index] = max(real(squeeze(theory_contrast_control.yy_freq(1,1,pos_indices,ii)))/norm_factor);
    max_freq(ii) = theory_contrast_control.params.omega(pos_indices(index))*1e3;
end

subplot(2,3,2); hold on;
indices_to_plot = [1:4];
tran_level = [1/4 1/4 1 1];
for kk = 1:length(indices_to_plot)
    ii = indices_to_plot(kk);
    norm_factor = real(squeeze(theory_contrast_light.yy_freq(1,1,theory_contrast_light.params.ind0,ii)));
    if ii == 4 || ii == 2
        gg(ii) = plot(theory_contrast_light.params.omega*1e3,(real(squeeze(theory_contrast_light.yy_freq(1,1,:,ii)))/norm_factor),...
            '-','linewidth',1.5,'color',[1 0.0353 0 tran_level(ii)]);
        xlim([0 100])
    end
    
    pos_indices = find(theory_contrast_light.params.omega*1e3>0);
    [gamma_power_light(ii), index] = max(real(squeeze(theory_contrast_light.yy_freq(1,1,pos_indices,ii)))/norm_factor);
    max_freq(ii) = theory_size_light.params.omega(pos_indices(index))*1e3;
end
set(gca,'fontsize',16)
xlabel('Freq (Hz)')
ylabel('c_{e1e1}(f) (normalized)')
legend([hh(2) gg(2) hh(4) gg(4)],{'low','low (light)','high','high (light)'})

subplot(2,3,5); hold on
plot(gamma_power_control,'.-','markersize',16,'color','k',...
        'linewidth',1.5)
plot(gamma_power_light,'.-','markersize',16,'color',[1 0.0353 0],...
        'linewidth',1.5)
ylabel('Normalized gamma power')
xlabel('Contrast')
set(gca,'fontsize',16)
xticks([indices_to_plot(1):indices_to_plot(end)])
xticklabels({'low','','','high'})
xlim([indices_to_plot(1) indices_to_plot(end)])
ylim_temp = ylim;
ylim([1 ylim_temp(2)])

figure(32); 
subplot(1,3,2); hold on;
plot(gamma_power_control,gamma_power_light,'k.-','linewidth',2,'markersize',30)
x_temp = min(gamma_power_control):0.001:max(gamma_power_control);
plot(x_temp,x_temp,'k')
set(gca,'fontsize',16)
xlabel('Gamma Control')
ylabel('Gamma Light')
xlim([min(gamma_power_control) max(gamma_power_control)])
xticks([])
yticks([])



%%
clearvars gamma_power_control gamma_power_light hh gg
figure(11);
subplot(2,3,3); hold on;
indices_to_plot = [1:4];
for kk = 1:length(indices_to_plot)
    ii = indices_to_plot(kk);
    
    norm_factor = real(squeeze(theory_surround_control.yy_freq(1,1,theory_surround_control.params.ind0,ii)));
    
    if kk == 1 || kk == 4
        hh(ii)  = plot(theory_surround_control.params.omega*1e3,(real(squeeze(theory_surround_control.yy_freq(1,1,:,ii)))/norm_factor),...
        '-','linewidth',1.5,'color',[0 0 0 ii/length(indices_to_plot)]);
        xlim([0 100])
    end
    
    pos_indices = find(theory_surround_control.params.omega*1e3>0);
    [gamma_power_control(ii), index] = max(real(squeeze(theory_surround_control.yy_freq(1,1,pos_indices,ii)))/norm_factor);
    max_freq(ii) = theory_surround_control.params.omega(pos_indices(index))*1e3;
end
set(gca,'fontsize',16)
xlabel('Freq (Hz)')
ylabel('c_{e1e1}(f) (normalized)')
legend('small','large')

subplot(2,3,3); hold on;
indices_to_plot = [1:4];
for kk = 1:length(indices_to_plot)
    ii = indices_to_plot(kk);
    norm_factor = real(squeeze(theory_surround_light.yy_freq(1,1,theory_surround_light.params.ind0,ii)));
    if ii == 4 || ii == 1
        gg(ii) = plot(theory_surround_light.params.omega*1e3,(real(squeeze(theory_surround_light.yy_freq(1,1,:,ii)))/norm_factor),...
            '-','linewidth',1.5,'color',[1 0.0353 0 ii/4]);
        xlim([0 100])
    end
    
    pos_indices = find(theory_surround_light.params.omega*1e3>0);
    [gamma_power_light(ii), index] = max(real(squeeze(theory_surround_light.yy_freq(1,1,pos_indices,ii)))/norm_factor);
    max_freq(ii) = theory_surround_light.params.omega(pos_indices(index))*1e3;
end
set(gca,'fontsize',16)
xlabel('Freq (Hz)')
ylabel('c_{e1e1}(f) (normalized)')
legend([hh(1) gg(1) hh(end) gg(end)],{'cross','cross (light)','iso','iso (light)'})

subplot(2,3,6); hold on
max_gamma = max(gamma_power_control);
plot(flip(gamma_power_control),'.-','markersize',16,'color','k',...
        'linewidth',1.5)
plot(flip(gamma_power_light),'.-','markersize',16,'color',[1 0.0353 0],...
        'linewidth',1.5)
ylabel('Normalized gamma power')
xlabel('Surround')
set(gca,'fontsize',16)
xticks([indices_to_plot(1):indices_to_plot(end)])
xticklabels({'iso','','','cross'})
xlim([indices_to_plot(1) indices_to_plot(end)])
ylim_temp = ylim;
ylim([1.5 ylim_temp(2)])

figure(32); 
subplot(1,3,3); hold on;
plot(flip(gamma_power_control),flip(gamma_power_light),'k.-','linewidth',2,'markersize',30)
x_temp = min(gamma_power_control):0.001:max(gamma_power_control);
plot(x_temp,x_temp,'k')
set(gca,'fontsize',16)
xlabel('Gamma Control')
ylabel('Gamma Light')
xlim([min(gamma_power_control) max(gamma_power_control)])
xticks([])
yticks([])

