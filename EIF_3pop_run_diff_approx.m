%%
% This runs the spiking network simulation using mex code.
%
% Saves the full spike train of all trials, and then condenses the data to
% find the average steady state firing rate, the power-spectrum and the
% cross-spectrum
%
% Denote the desired location to save the large data file (i.e., all spikes)
% in data_dir
%
% Meant for small simulations to simply check the steady state firing rates
% against the mean-field theory
%
% Written by Gregory Handy, 08/11/2021
%%
clear; clc; close all;

restoredefaultpath;
folder = fileparts(which('EIF_3pop_run_diff_approx.m')); 
addpath(genpath(folder));
rmpath(folder)
rmpath('./Cluster_code')


%% Compile the mex code

% Use to compile the mex code (do after each edit to the mex file)
% Run the following if running mex for the first time
% run mex -setup -c
mex ./Mex_Functions/EIF_mex_diff_approx_delay.c

%% Save the data in this directory
if ismac
    data_dir = '/Users/gregoryhandy/Research_Local/EIF_Gamma_Data/';
else
    data_dir = '/user_data/ghandy/EIF_Gamma_Data/';
end

%% Load the parameters
light = {'range',0};
surround_diff = 0;
params = EIF_params_official_fn(surround_diff,1,light,'med');

%% Seed the random number generator
rng('shuffle');
s = rng();
random_seed = s.Seed;
% still random, but saves the spikes in this designated spot
% random_seed = 1002; 

%% Create the connectivity matrix to be used for all simulations
tic;
fprintf('Creating the connectivity matrix \n');
[wind,wipost,wstr] = ...
    Diff_approx_gen_weights(params.Ncells,params.p,params.J);

% offset the index by one for the MEX code (C starts indexing at 0)
wind_mex = int32(wind-1);
wipost_mex = int32(wipost-1);
pinds_mex = int32(params.pinds-1);
toc;

%% Conditions to loop over
r_vip = [1:(12-1)/7:12];
num_conds = length(r_vip);

%%
tic;
fprintf('Running the spiking simulations \n');
parfor cond_num = 1:num_conds
     
    %% Load the parameters
    params = EIF_params_official_fn(surround_diff,1,{'range',r_vip(cond_num)},'med'); %#ok<PFTUSW>
   
    %% Define the feedforward input (mean and variance)
    mu_ext = params.mu_bg + params.mu_stim + params.mu_vip;
    sigma_ext = sqrt(params.sigma_bg.^2+params.sigma_stim.^2+...
        params.sigma_vip.^2).*sqrt(2*params.tau_m);
         
    %% Run the simulation
    [rates_trial,times,tinds] = ...
        EIF_mex_diff_approx_delay(params.T, params.NT, params.recstart,...
        params.Npop, params.Ntot, params.Ncells, params.dt, params.rates, params.p, params.J,...
        params.tau_s, params.tau_m, params.EL, params.vth, params.vre, params.tau_ref,...
        wind_mex, wipost_mex, wstr, pinds_mex,random_seed,mu_ext,...
        sigma_ext,params.DeltaT,params.VT, params.tauDelay,...
        params.maxDelay,params.shared_noise_mag.*sqrt(2*params.tau_m),...
        params.maxns);   

    save_data(data_dir,random_seed,rates_trial,times,tinds,params,cond_num)
end
toc;

%% Find the steady state firing rates
tic;
rates_ave = zeros(params.Npop,num_conds);
for ii = 1:num_conds
    
    %% Load the data
    fprintf('Load data from stim %d \n',ii)
    
    file_name = sprintf('EIF_stim_num_%d_%d',ii,random_seed);
    name_full = strcat(data_dir,file_name);
    load(name_full,'rates_trial','times','tinds')
    
    %% Estimate the firing rate
    fprintf('Estimating firing rates \n')
    for jj = 1:params.Npop
        rates_ave(jj,ii) = mean(rates_trial(1,params.pinds(jj):params.pinds(jj+1)-1));
    end
    
end
toc;

%% Plot against the mean-field theory
color_scheme = [0.8500, 0.3250, 0.0980; 0, 0.4470, 0.7410; 0.4660, 0.6740, 0.1880;...
    0.8500, 0.3250, 0.0980; 0, 0.4470, 0.7410; 0.4660, 0.6740, 0.1880]; 

mf_data = load(sprintf('./Data_sets/EIF_LR_MF_data_vip_surround_diff_%d',surround_diff));

figure(1); clf; hold on
for ii = 1:6
    plot(mf_data.rates_trial_ave(ii,:),'linewidth',1.5,'color',color_scheme(ii,:))
    if ii <=3
        plot(rates_ave(ii,:),'.','markersize',16,'color',color_scheme(ii,:))
    else
        plot(rates_ave(ii,:),'o','markersize',10,'color',color_scheme(ii,:))
    end
end



