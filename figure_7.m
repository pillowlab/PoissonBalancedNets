%% Figure 7
%
% =======
% Compares the local and population framework for two simulations, a 1-D perfect
% integrator and a 2-D harmonic oscillator. The first part of the code 
% outputs the raster plots and the read-outs for all four conditions. The 
% second part outputs plots comparing the read-out error and spikes fired 
% for the four conditions.
% =======
%
% Dependencies: local_framework.m, population_framework.m, test_approx.m
%
%% Set up stimulus
NT = 1; % max time
dt = .1/1000; % dt 
tt = (dt:dt:NT); % time vector
nt = length(tt); % number of time bins

%create stimulus
blocksize = .1; % size of stimulus blocks 
nblock = round((blocksize/dt)); %number of stimulus blocks
stimrnge = [-10 10]; %range of stimulus intensities 
stim = reshape(ones(nblock,1)*rand(1,nt/nblock),[],1)'; % generate random stimulus
stim = stim*diff(stimrnge)+stimrnge(1); % shift to have desired range
stim(1:nblock/10) = 0; % set first block to zeros

%plot stimulus
figure(1)
subplot(521);
h = plot(tt, tt*0, 'k--', tt, stim);
set(h(2), 'linewidth', 2);
xlabel('Time (s)');
set(gca, 'FontSize',15, 'linewidth',1)
ylabel('stimulus')
box off; drawnow;

%% Set up the population of GLM neurons
params.N = 200; % number of neurons
params.wmean = .02; % mean weight value 
params.wsig = .01; % stdev of noise in weights
params.taud = .2; % rate decay for filtered spike train
xtau = 0.02; % time constant for x dynamics
params.A = -1/xtau; % x dynamics
params.mu = 0; % quadratic cost on spiking
%local framework parameters
params.fmax = 100; % max firing rate
params.alpha = 800; % precision 
params.fmin = 0; % baseline firing rate
%population framework parameters
params.beta = 25; % scaling for pspike

%% Run simulations

% local
params.tdel = 0; % time delay in bins
[sps,xh,xx] = local_framework(stim, dt, params); %run local framework simulation
params.t_del = 50; %time delay in bins 
[sps2,xh2,~] = local_framework(stim, dt, params); %run local framework simulation

% population
params.t_del = 0;% time delay in bins
[sps3,xh3,~] = population_framework(stim, dt, params); %run population framework simulation
params.t_del = 50; % time delay in bins
[sps4,xh4,~] = population_framework(stim, dt, params); %run population framework simulation

%% Plot spikes and read-out
% raster plot for local without delay
subplot(5,2,[3])
N = params.N; % number of neurons
[iiinh,jjinh] = find(sps(1:N/2,:)); % spikes of inhibitory
[iiexc,jjexc] = find(sps(N/2+1:N,:)); % spikes of excitatory
h = plot(jjinh*dt,iiinh, '.', jjexc*dt,iiexc+N/2,'.');
ylabel('Neuron')
set(gca,'xlim', [0 NT],'ylim',[0 N]); box off
set(gca, 'FontSize',15, 'linewidth',1)
title('local');

%rarster plot for local with delay
subplot(5,2,[5])
[iiinh,jjinh] = find(sps2(1:N/2,:)); 
[iiexc,jjexc] = find(sps2(N/2+1:N,:)); 
h = plot(jjinh*dt,iiinh, '.', jjexc*dt,iiexc+N/2,'.');
set(gca,'xlim', [0 NT],'ylim',[0 N]); box off
set(gca, 'FontSize',15, 'linewidth',1)

%raster plot for population without delay
subplot(5,2,[4])
[iiinh,jjinh] = find(sps3(1:N/2,:));
[iiexc,jjexc] = find(sps3(N/2+1:N,:));
h = plot(jjinh*dt,iiinh, '.', jjexc*dt,iiexc+N/2,'.');
set(gca,'xlim', [0 NT],'ylim',[0 N]); box off
set(gca, 'FontSize',15, 'linewidth',1)
title('population');

%raster plot for population with delay
subplot(5,2,[6])
[iiinh,jjinh] = find(sps4(1:N/2,:)); 
[iiexc,jjexc] = find(sps4(N/2+1:N,:)); 
h = plot(jjinh*dt,iiinh, '.', jjexc*dt,iiexc+N/2,'.');
set(gca,'xlim', [0 NT],'ylim',[0 N]); box off
set(gca,'xlim', [0 NT],'ylim',[0 N])
box off
set(gca, 'FontSize',15, 'linewidth',1)

% read-out for local with and withou delay
subplot(5,2,[7 9])
h = plot(tt, xx, 'k-', tt, xh, tt, xh2,'-',tt, tt*0,'k--');
set(h(1),'linewidth',1);
legend('target', 'no delay', 'delay','location', 'northwest')
legend boxoff 
box off
set(gca, 'FontSize',15, 'linewidth',1)
xlabel('Time (s)');
set(gca,'ylim',[min(xx)+.4*min(xx) max(xx)+.4*max(xx)]);

% read-out for population with and without delay
subplot(5,2,[8 10])
h = plot(tt, xx, 'k-', tt,xh3, tt,xh4,'-',tt, tt*0,'k--');
set(h(1),'linewidth',1);
legend('targ', 'no delay', 'delay', 'location', 'northwest')
legend boxoff 
box off
set(gca, 'FontSize',15, 'linewidth',1)
set(gca,'ylim',[min(xx)+.4*min(xx) max(xx)+.4*max(xx)]);
set(gcf, 'renderer', 'Painters')

%% Calculate errors in local and population frameworks due to time delay
% runs a simulation for a given time delay for each framework a certain 
% number of times (runs) and gets the average. compares that to a minimum
% error in encoding x that arises from the mismatch of using a time-delayed
% euler approximation (using test_approx.m)

runs = 30; %number of times to repeat the local and population simulations
tdels = 0:5:45; % array of time delays to iterate over 

%initialize variables
err_min = zeros(1,length(tdels)); %mimimum error from time delay
err_net = zeros(1,length(tdels));%mean error in local framework over all runs
err_loc = zeros(1,length(tdels));%error in population framework over all runs

for i = 1:length(tdels)
    
params.tdel = tdels(i); % get value for time delay
[xx,xt] = test_approx(stim,dt,params); %run simulation to get theoretical bound on error
err_min(i) = 1 - sum((xx-xt).^2)/sum((xx - mean(xx)).^2); %bound on R^2 value for this time delay

loc = zeros(1,runs); % R^2 value for all of the runs, local framework
net = zeros(1,runs);% R^2 value for all of the runs, population framework

for j = 1:runs
[~,xh1,x1] = local_framework(stim, dt, params); % run local framework
[~,xh2,x2] = population_framework(stim, dt, params); % run population framework

loc(j) = 1 - sum((x1-xh1).^2)/sum((xh1 - mean(xh1)).^2); %error in local framework read-out 
net(j) = 1 - sum((x2-xh2).^2)/sum((xh2 - mean(xh2)).^2); %error in population framework read-out

end

err_loc(i) = mean(loc); %get average error for all local framework simulations with a given time delay
err_net(i) = mean(net); %get average error for all population framework simulations with a given time delay

end
%% Plot errors 
figure(2)
plot(tdels, err_min) % minimum error over all time delays
hold on
plot(tdels, err_net) % error in population framework over all time delays
hold on
plot(tdels, err_loc) % error in local framework over all time delays
set(gca, 'ylim', [0.4 max(err_min(1:7))], 'xlim', [0 45], 'linewidth', 1, 'FontSize', 15)
box off
legend('min', 'global', 'local')
legend boxoff
xlabel('time delay (bins)')
ylabel('R^2 value')
title('decoding error')
