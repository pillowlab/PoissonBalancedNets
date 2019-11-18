%% Figure 8
%
% =======
% Runs a long simulation of a randum stimulus in order to compute across
% and within population cross-correlations and autocorrelations for the
% BSN model and the population and local frameworks. 
% =======
% 
% Dependencies: local_framework.m, population_framework.m, bsn.m
% 
%% set up stimulus
%load thatsgood
NT = 80; % max time (s)
dtStim = .0001; % dt in s
time = (dtStim:dtStim:NT); % time vector
nt = length(time); %number of time steps

% generate random stimulus
stim = 50*randn(1, nt); %white noise
fc = 500; %cutoff frequency
stim2 = lowpass(stim,fc, 1/dtStim); %low pass filter

%plot stimulus
figure(1)
subplot(211);
plot(time,stim2)
title('low pass filtered white noise stimulus')

subplot(212);
plot(time, time*0, 'k--', time,cumsum(stim2))
%% Set up the population of GLM neurons
N = 40;
params.N = N; % Number of neurons
params.taud = 0.5; %decay of filtered spike train 
params.mu = 0; % quadratic cost on spiking
params.tdel = 5; % time delay (in bins)
params.wsig = .01; % noise in weights 
params.A = 1/Inf; % integrator
params.wmean = 1/N*ones(1,size(params.A,1))/10; % mean weight value

%% Run local framework
% NOTE: The local framework typically does not generate enough spikes to
% get an accurate estimate of the cross- and auto-correlations. To increase
% the number of spikes, increase f_max or alpha. 

%local framework parameters
params.fmax = 5000; % max firing rate
params.alpha = 2000; %  precision (higher -> more deterministic spiking)
params.fmin = 5; % 2 baseline log firing rate (spikes/bin)

%run local framework
[o, ~, ~] = local_framework(stim2, dtStim, params);
oLocal = o;

%% population framework
%population framework parameters
params.beta = 20; % scaling of pspike for population simulation

% run population simulation
[o, ~, ~] = population_framework(stim2, dtStim, params);
oPopulation = o;

%% BOERLIN MODEL
%NOTE: if these parameter settings are not giving enough spikes,
%lower the cost on spiking (params.nu and .mu)  

%BSN parameters
params.dt = dtStim; %dt
params.sp_cond = 1; %spiking condition (1 = randomly select one to fire, 0 is let all fire)
params.tau_v = .2; %voltage decay
params.nu = 0; %linear cost orig 
params.mu = 0; %quadratic cost 
params.sigv =0; %voltage noise

%run BSN simulation
[o, ~, ~] = bsn(stim,params);
oBSN= o;
%%  Compute and plot cross-correlations

% TOGGLE THIS TO GET XCORR AND AUTOCORR FOR DIFFERENT MODELS
%sps = oLocal';
sps = oPopulation';
%sps = oBSN';

meanSC = mean(sum(sps)); % mean spike count
fprintf('mean spike count = %.1f\n', meanSC);

maxlag = 50;
xce = zeros(maxlag*2+1,1); % cross-correlation excitatory pop
xci = zeros(maxlag*2+1,1); % cross-correlation inhibitory pop
xcx = zeros(maxlag*2+1,1); % cross-correlation across pop
ace = zeros(maxlag*2+1,1); % auto-correlation excitatory pop
aci = zeros(maxlag*2+1,1); % auto-correlation inhibitory pop

% compute cross- and auto-correlations
for ii = 1:N/2
    ace = ace + xcorr(sps(:,ii), sps(:,ii), maxlag, 'unbiased');
    aci = aci + xcorr(sps(:,ii+N/2), sps(:,ii+N/2), maxlag, 'unbiased');        
    for jj = ii+1:N/2
    xce = xce + xcorr(sps(:,ii), sps(:,jj), maxlag, 'unbiased');
    xci = xci + xcorr(sps(:,ii+N/2), sps(:,jj+N/2), maxlag, 'unbiased');
    xcx = xcx + xcorr(sps(:,ii), sps(:,jj+N/2), maxlag, 'unbiased');
    end
end

%set values at origin to NaN
ace(maxlag+1) = NaN;
aci(maxlag+1) = NaN;

%% plot auto and cross correlation for one of the frameworks
tt = (-maxlag:maxlag)*dtStim*1000; %for plotting
figure(1)
clf;
% autocorrellation for both populations
subplot(311)
plot(tt, (ace+aci)/2, 'linewidth', 1.5)%-mean(B.ace+B.aci)], 'linewidth', 1.5)
ylabel('Auto')
box off
set(gca, 'FontSize', 14, 'linewidth', 1)

%cross correlation within populatoins
subplot(312)
plot(tt, xce, 'linewidth', 1.5)
hold on
plot(tt, xci, 'linewidth', 1.5)
ylabel('within')
set(gca, 'FontSize', 14, 'linewidth', 1)
box off
legend('excitatory', 'inhibitory')
legend boxoff

%cross correlation across populations
subplot(313)
plot(tt, xcx, 'linewidth', 1.5)
set(gca, 'FontSize', 14, 'linewidth', 1)
box off
ylabel('across')
xlabel('Time lag (ms)');