%% Figure 2 
%
% =======
% Replicates figure 1.C from Boerlin et al. (2013), which shows the
% BSN integrating a noisy stimulus. This code outputs the read-out for
% two different spiking condtions: allowing only one neuron over threshold 
% to fire or allowing all over threshold to fire.
% =======
%
% Dependencies: bsn.m
% 
%% Set up stimulus
tmax = 1.5; % length of simulation 
params.dt = .1/1000; % time bin
t = (params.dt:params.dt:tmax); % time vector

stim = zeros(1,length(t)); % stimulus
stim(t>.25&t<.55) = 50;
stim(t>.8&t<1) = -100;
stim = stim + 0.01*randn(1,length(t));

%% Set up simulation parameters
params.N = 400; %number of neurons
params.wmean = 0.1; % mean decoding weight
params.wsig = 0; % weight variance
params.nu = 10e-5; %linear cost orig 
params.mu = 10e-6; %10e-4 %quadratic cost 
params.sigv = 10e-4; %voltage noise
params.taud = 1/10; % sp train decay time constant
params.tau_v = 1/20; % voltage decay time constant

% Set spiking condition
% 0 = all above threshold fire
% 1 = randomly select from neurons above thresh
params.sp_cond = 0;
%% Run simulation
[ss, xh, xx] = bsn(stim,params);

%% Plot results
figure(1)
% Plot stimulus
subplot(5, 2, 1)
plot(t, stim,'linewidth', 2)
box off
set(gca, 'FontSize',15, 'linewidth',1)
xlabel('Time (s)')
ylabel('Stimulus intensity')

% Raster plot 
a = subplot(5, 2, [3 5])
cla(a)
[iiinh,jjinh] = find(ss(1:params.N/2,:)); % spikes of inhibitory
[iiexc,jjexc] = find(ss(params.N/2+1:params.N,:)); % spikes of excitatory
plot(jjinh*params.dt,iiinh, '.', jjexc*params.dt,iiexc+params.N/2,'.');
set(gca, 'ylim',[1 params.N]);
box off
set(gca, 'FontSize',15, 'linewidth',1)
ylabel('Neurons')
title('All fire')
%Read out and estimate
a = subplot(5, 2, [7 9])
cla(a)
plot(t,xx, 'k','linewidth', 1.25)
hold on
plot(t,xh, '--','linewidth', 2);
set(gca, 'ylim', [-10 20], 'xlim', [0 1.5], 'FontSize', 12)
box off
set(gca, 'FontSize',15, 'linewidth',1)
ylabel('Read-out')
%% Randomly select neuron to fire
params.sp_cond = 1; % 0 let all fire, 1 random select firing, 2 max over T fires
[ss, xh, xx] = bsn(stim,params);

%% Plot
% Raster plot
a = subplot(5, 2, [4 6])
cla(a)
[iiinh,jjinh] = find(ss(1:params.N/2,:)); % spikes of inhibitory
[iiexc,jjexc] = find(ss(params.N/2+1:params.N,:)); % spikes of excitatory
plot(jjinh*params.dt,iiinh, '.', jjexc*params.dt,iiexc+params.N/2,'.');
set(gca, 'ylim',[1 params.N]);
box off
set(gca, 'FontSize',15, 'linewidth',1)
title('Select one')

%Read out and target
a = subplot(5, 2, [8 10])
cla(a)
plot(t,xx, 'k','linewidth', 1.25)
hold on
plot(t,xh, '--','linewidth', 2);
set(gca, 'ylim', [-10 20], 'xlim', [0 1.5], 'FontSize', 12)
box off
set(gca, 'FontSize',15, 'linewidth',1)
