clear; close all; clc;
%%
restoredefaultpath;
folder = fileparts(which('plot_paper_raster.m')); 
addpath(genpath(folder));
rmpath(folder)

%%
data_dir = '/Users/gregoryhandy/Research_Local/EIF_Gamma_Data/';
% random_seed = 1377525081; % cross (small dt)
random_seed = 1635400248;
stim_num = 2;

file_name = sprintf('EIF_stim_num_%d_%d',stim_num,random_seed);


load(strcat(data_dir,file_name))

%%

figure(1); clf;
subplot(3,1,2:3)
plot_raster_v2(times, tinds, params.Ntot, 3, params.Ncells, params.pinds, 1 )

subplot(3,1,2:3)
tStart = 893;
tEnd = 893.25;
xlim([tStart tEnd])
xticks([tStart:0.25:tEnd])
xticklabels([0 250 500])
xlabel('Time (msec)')
ylabel('Neuron Index')
box off
%%
clearvars temp;
subplot(3,1,1)
delta_t = 0.005;
tvec = [tStart:delta_t:tEnd];

times_small = times(find(times/1000>tStart & times/1000<=tEnd & tinds < params.Ncells(1)));
tinds_small = tinds(find(times/1000>tStart & times/1000<=tEnd & tinds < params.Ncells(1)));

temp = zeros(length(tvec)-1,1);
tic;
for i = 1:length(tvec)-1
    temp(i) = length(find(times_small/1000>tStart+delta_t*(i-1) & times_small/1000<=tStart+delta_t*i & tinds_small < params.Ncells(1)));
end
toc;
plot(tvec(1:end-1),temp/(params.Ncells(1)*delta_t),'color',[38 30 101]/255,'linewidth',1.5)
box off
xlim([tStart tEnd])
xticks([])
% xticklabels([0 500 1000])
ylabel('Firing rate (Hz)')
set(gca,'fontsize',16)

mean(rates_trial(1:params.Ncells(1)))



