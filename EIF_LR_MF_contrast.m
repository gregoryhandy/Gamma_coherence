%%
% Finds the cross-correlation function for multiple freq. and time lags
% using mean field theory for a range of contrast values 
%
% Written by Gregory Handy, 08/24/2021
%%
clear; close all; clc;

restoredefaultpath;
folder = fileparts(which('EIF_LR_MF_contrast.m')); 
addpath(genpath(folder));
rmpath(folder)
tic;
%% Compile mex code (if necessary)

mex ./Mex_functions/calc_Rate_cuEIF.cpp
mex ./Mex_functions/calc_Power_cuEIF.cpp
mex ./Mex_functions/calc_Susc_cuEIF.cpp

%% Load parameters

light = {'off'}; stim_size = 'med';
params = EIF_params_official_fn(0,1,light,stim_size);

%% Preallocate 
contrast = [0.5:(1-0.5)/3: 1];
num_conds = length(contrast);

rates_trial = zeros(num_conds,params.Npop);
yy_freq = zeros(params.Npop,params.Npop,params.bins,num_conds);
yy_time = zeros(params.Npop,params.Npop,params.bins,num_conds);

%% Loop over the stimuli
parfor jj = 1:num_conds
    
    % Load the parameters
    params = EIF_params_official_fn(0,contrast(jj),light,stim_size); %#ok<PFTUSW>
    
    % Compute the linear response theory
    [rates_trial(jj,:),yy_freq(:,:,:,jj),yy_time(:,:,:,jj)]=EIF_linear_response_fn(params)
end
rates_trial_ave = rates_trial'*10^3; % convert the rates to Hz

%% Save the results
% filename = sprintf('./Data_sets/Exp_compare/EIF_LR_MF_data_contrast_light_%s',light{1});
% save(filename,'rates_trial_ave','params','yy_time','yy_freq')

%%
color_scheme = [0.8500, 0.3250, 0.0980; 0, 0.4470, 0.7410; 0.4660, 0.6740, 0.1880;...
    0.8500, 0.3250, 0.0980; 0, 0.4470, 0.7410; 0.4660, 0.6740, 0.1880]; 

figure(1); clf;
for ii = 1:params.Npop
    h(ii) = plot(rates_trial_ave(ii,:),'.-','markersize',16,'linewidth',1.5,'color',color_scheme(ii,:));
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
indices_to_plot = [1:num_conds];
for kk = 1:length(indices_to_plot)
    ii = indices_to_plot(kk);
    norm_factor = 1+0*real(squeeze(yy_freq(1,1,params.ind0,ii)));
    plot(params.omega*1e3,(real(squeeze(yy_freq(1,1,:,ii)))/norm_factor),...
        '-','linewidth',1.5,'color',[0 0 0 ii/num_conds])
    xlim([0 100])
    
    pos_indices = find(params.omega*1e3>0);
    [val, index] = max(real(squeeze(yy_freq(1,1,pos_indices,ii))));
    max_freq(ii) = params.omega(pos_indices(index))*1e3;
end
set(gca,'fontsize',16)
xlabel('Freq (Hz)')
ylabel('c_{e1e1}(f)')


%%
clearvars h;
figure(2); clf; hold on;
subplot(1,2,1); hold on

for kk = 1:length(indices_to_plot)
    ii = indices_to_plot(kk);
    norm_factor = real(squeeze(yy_freq(1,1,params.ind0,ii)));
    h(kk) = plot(params.omega*1e3,(real(squeeze(yy_freq(1,1,:,ii)))/norm_factor),...
        '-','linewidth',1.5,'color',[0 0 0 ii/num_conds]);
    xlim([0 100])
    
    pos_indices = find(params.omega*1e3>10);
    [gamma_power(ii), index] = max(real(squeeze(yy_freq(1,1,pos_indices,ii)))/norm_factor);
    max_freq(ii) = params.omega(pos_indices(index))*1e3;
    
end
set(gca,'fontsize',16)
xlabel('Freq (Hz)')
ylabel('Power (normalized)')
legend([h(1),h(end)],{'Cross','Iso'}) 

subplot(1,2,2); hold on
plot(gamma_power,'.-','linewidth',1.5,'markersize',16,'color','k')
ylabel('Normalized gamma power')
xlabel('Contrast')
set(gca,'fontsize',16)
xticks([indices_to_plot(1) indices_to_plot(end)])
xticklabels({'Cross','Iso'})
xlim([indices_to_plot(1) indices_to_plot(end)])



