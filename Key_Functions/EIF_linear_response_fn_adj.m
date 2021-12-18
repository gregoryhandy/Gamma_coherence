function [rates,yy_freq,yy_time,rates_min] = EIF_linear_response_fn_adj(params)

% Create the connectivity matrix
C = params.p'.*params.Ncells'.*params.W;

%% Solve for the rates by fixed point iteration
E0 = params.EL +params.mu_bg + params.mu_stim+params.mu_vip;
sigma = sqrt(params.sigma_bg.^2+params.sigma_stim.^2+params.sigma_vip.^2 ...
    +params.shared_noise_mag.^2);

% Generate an initial estimate for the FP iteration
r0 = zeros(params.Npop,1);
rates = zeros(params.Npop,1);
for pc = 1:params.Npop
    r0(pc) = calc_Rate_cuEIF(E0(pc),sigma(pc),params.tau_ref(pc),...
        params.vre(pc),params.VT(pc),params.vth(pc),params.DeltaT(pc),...
        params.tau_m(pc),params.vlb,params.dv);
    
    rates(pc,1) = r0(pc);
end
rates_temp = zeros(params.Npop,1);

% Fixed point iteration
iter_num = 1;
max_change = 10;
rates_saved = zeros(10000-1,6);
while max_change > 1e-7 && iter_num < 10000
    for pc = 1:params.Npop
        E0_eff = E0(pc)+C(pc,:)*rates;
        sigma_eff = sqrt(sigma(pc).^2);
        
        rates_temp(pc) = calc_Rate_cuEIF(E0_eff,sigma_eff,...
            params.tau_ref(pc),params.vre(pc),params.VT(pc),params.vth(pc),...
            params.DeltaT(pc),params.tau_m(pc),params.vlb,params.dv);
    end
    
    max_change = max(abs(rates-rates_temp));
    rates = rates_temp;
    rates_saved(iter_num,:) = rates_temp;
    iter_num = iter_num + 1;
end


if iter_num == 10000
    iter_num
    rates_saved = rates_saved(9000:10000-1,:);
    
    [~,max_index] = max(rates_saved(:,1));
    rates = rates_saved(max_index,:)';
    
    [~,min_index] = min(rates_saved(:,1));
    rates_min = rates_saved(min_index,:)';
else
    rates_min=nan(6,1);
end

%% Calculate the uncoupled power spectrum, susceptibility, and Fourier
% transformed synaptic kernel for every cell in the network at the
% frequency values w.
Ct0 = zeros(params.Npop,params.bins);
At = zeros(params.Npop,params.bins);
Ft = zeros(params.Npop,params.bins);

yy = zeros(params.Npop,params.Npop,params.bins);
K = zeros(params.Npop,params.Npop);
yy0 = zeros(params.Npop,params.Npop);
I = eye(params.Npop);

for pc = 1:params.Npop
    % Calculate the effective rest potential in the network
    E0_eff = E0(pc) + C(pc,:)*rates;
    sigma_eff = sigma(pc);
    
    Ct0(pc,:) = calc_Power_cuEIF(params.omega,E0_eff,sigma(pc),params.tau_ref(pc),...
        params.vre(pc),params.VT(pc),params.vth(pc),params.DeltaT(pc),...
        params.tau_m(pc),params.vlb,params.dv,rates(pc));
    
    At(pc,:) = calc_Susc_cuEIF(params.omega,E0_eff,sigma(pc),params.tau_ref(pc),...
        params.vre(pc),params.VT(pc),params.vth(pc),params.DeltaT(pc),...
        params.tau_m(pc),params.vlb,params.dv,rates(pc));
    
    % Calculate the values of Fourier transform of the synaptic kernel for cell i.
    for j = 1:params.bins
        Ft(pc,j) = exp(1i*-2*pi*params.omega(j)*params.tau_d(pc))/((1+1i*2*pi*params.omega(j)*params.tau_s(pc))^2);
    end
end

%% Solve the linear response equations for the auto- and cross-spectra of
% every pair of cells in the network and store in yy.
for j = 1:params.bins
    
    for k = 1:params.Npop
        for l = 1:params.Npop
            K(k,l) = C(k,l)*At(k,j)*Ft(l,j);
        end
        yy0(k,k) = Ct0(k,j);
    end
    
    yy0_ave = yy0./params.Ncells;

    global_noise_term = diag(At(:,j))*(sqrt(2*params.tau_m).*params.shared_noise_mag)'*...
        (diag(At(:,j))*(sqrt(2*params.tau_m).*params.shared_noise_mag)')';    
    global_correction_term = ((sqrt(2*params.tau_m).^2.*params.shared_noise_mag.^2).*...
        abs(diag(At(:,j))).^2)./params.Ncells;
    
    yy(:,:,j) = (I-K)\(yy0_ave+global_noise_term-global_correction_term)/(I-K');
end


%% Find the spike train cross-covariance function C(s) by inverse Fourier transformation
% Calculate the inverse Fourier transform of the auto-/cross-spectra for
% every pair.
yy_freq = zeros(params.Npop,params.Npop,params.bins);
yy_time = zeros(params.Npop,params.Npop,params.bins);
for k = 1:params.Npop
    for l = 1:params.Npop
        temp = squeeze(yy(k,l,:));
        yy_freq(k,l,:) = temp;
        
        [~,temp_ccg] = inv_f_trans_on_vector(params.omega,temp);
        yy_time(k,l,:) = real(temp_ccg);
    end
end

end

