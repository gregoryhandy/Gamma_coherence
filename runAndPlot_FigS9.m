%%
% Finds the cross-correlation function for multiple freq. and time lags
% using mean field theory for a range of contrast values for the modified
% model where VIP is modulated by the stimulus
%
% Corresponds to Fig. S9 in:
%   J Veit, G Handy, DP Mossing, B Doiron, H Adesnik. 
%   Cortical VIP neurons locally control the gain but globally control 
%   the coherence of gamma band rhythms.
%
% Written by Gregory Handy, 08/24/2021
%%
clear; close all; clc;

restoredefaultpath;
folder = fileparts(which('runAndPlot_FigS9.m')); 
addpath(genpath(folder));
rmpath(folder)
%% Compile mex code (if necessary)

% mex ./Mex_functions/calc_Rate_cuEIF.cpp
% mex ./Mex_functions/calc_Power_cuEIF.cpp
% mex ./Mex_functions/calc_Susc_cuEIF.cpp

%% Load parameters

light = {'off'}; stim_size = 'med';
params = EIF_params_official_fn(0,1,light,stim_size);

%% Preallocate 
contrast = [0.5:(1-0.5)/3: 1];
num_conds = length(contrast);

rates_trial = zeros(num_conds,params.Npop);
yy_freq = zeros(params.Npop,params.Npop,params.bins,num_conds);
yy_time = zeros(params.Npop,params.Npop,params.bins,num_conds);

rates_trial_r2 = zeros(num_conds,params.Npop);
yy_freq_r2 = zeros(params.Npop,params.Npop,params.bins,num_conds);
yy_time_r2 = zeros(params.Npop,params.Npop,params.bins,num_conds);

%% Loop over the stimuli
parfor jj = 1:num_conds
    
    % Load the parameters
    params = EIF_params_official_fn(0,contrast(jj),light,stim_size); %#ok<PFTUSW>
    
    % Compute the linear response theory
    [rates_trial(jj,:),yy_freq(:,:,:,jj),yy_time(:,:,:,jj)]=EIF_linear_response_fn(params)
end
rates_trial_ave = rates_trial'*10^3; % convert the rates to Hz

%%
parfor jj = 1:num_conds
    
    % Load the parameters
    params = EIF_params_official_fn_FigS9(0,contrast(jj),light,stim_size); %#ok<PFTUSW>
    
    % Compute the linear response theory
    [rates_trial_r2(jj,:),yy_freq_r2(:,:,:,jj),yy_time_r2(:,:,:,jj)]=EIF_linear_response_fn(params)
end
rates_trial_ave_r2 = rates_trial_r2'*10^3; % convert the rates to Hz

%%
color_scheme = [0, 0, .54; 0, .392, 0; 1, .749, 0;.545, 0, .545]; 
figure(1); clf; h = [];
subplot(1,3,1); hold on;
for ii = 1:3
    h(ii) = plot(rates_trial_ave(ii,:),'-','markersize',16,'linewidth',1.5,'color',color_scheme(ii,:));
    plot(rates_trial_ave_r2(ii,:),'--','markersize',16,'linewidth',1.5,'color',color_scheme(ii,:));
end

vipRates_r2 = [6 5 4 3];
h(4) = plot([4 4 4 4],'-','markersize',16,'linewidth',1.5,'color',color_scheme(4,:));
plot(vipRates_r2,'--','markersize',16,'linewidth',1.5,'color',color_scheme(4,:))

set(gca,'fontsize',16)
legend(h(1:4),{'E','PV','SST','VIP'})
ylabel('Firing rate (Hz)')
xticks([1 4])
xticklabels({'Low','High'})
xlabel('Contrast')

%%
indices_to_plot = [1:num_conds];
clearvars h;
subplot(1,3,2); hold on

for kk = 1:length(indices_to_plot)
    ii = indices_to_plot(kk);
    norm_factor = real(squeeze(yy_freq(1,1,params.ind0,ii)));
    h(kk) = plot(params.omega*1e3,(real(squeeze(yy_freq(1,1,:,ii)))/norm_factor),...
        '-','linewidth',1.5,'color',[0 0 0 ii/num_conds]);
    xlim([0 100])
    
    
    norm_factor_r2 = real(squeeze(yy_freq_r2(1,1,params.ind0,ii)));
    plot(params.omega*1e3,(real(squeeze(yy_freq_r2(1,1,:,ii)))/norm_factor_r2),...
        '--','linewidth',1.5,'color',[0 0 0 ii/num_conds]);
    
    pos_indices = find(params.omega*1e3>10);
    [gamma_power(ii), index] = max(real(squeeze(yy_freq(1,1,pos_indices,ii)))/norm_factor);
    max_freq(ii) = params.omega(pos_indices(index))*1e3;
    
    
    pos_indices = find(params.omega*1e3>10);
    [gamma_power_r2(ii), index] = max(real(squeeze(yy_freq_r2(1,1,pos_indices,ii)))/norm_factor_r2);
    max_freq_r2(ii) = params.omega(pos_indices(index))*1e3;
    
end
set(gca,'fontsize',16)
xlabel('Freq (Hz)')
ylabel('Power (normalized)')
legend([h(1),h(end)],{'Low','High'}) 

subplot(1,3,3); hold on
plot(gamma_power,'-','linewidth',1.5,'markersize',16,'color','k')
plot(gamma_power_r2,'--','linewidth',1.5,'markersize',16,'color','k')
ylabel('Normalized gamma power')
xlabel('Contrast')
set(gca,'fontsize',16)
xticks([indices_to_plot(1) indices_to_plot(end)])
xticklabels({'Low','High'})
xlim([indices_to_plot(1) indices_to_plot(end)])
legend('Untuned VIP','Tuned VIP')



