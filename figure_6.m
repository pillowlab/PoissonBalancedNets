%% Figure 6
%
% =======
% Compares the local and population framework for two simulations, a 1-D perfect
% integrator and a 2-D harmonic oscillator. The first part of the code 
% outputs the raster plots and the read-outs for all four conditions. The 
% second parrt outputs plots comparing the read-out error and spikes fired 
% for the four conditions.
% =======
% 
% Dependencies: local_framework.m, population_framework.m
% 
%% Set up 2D random stimulus
NT = 3; % simulation length
dtStim = .1/1000; % dt 
time = (dtStim:dtStim:NT); %time vector
nt = length(time); %number of time bins

% set up a stimulus consisting of random blocks of set length within a 
% certain y value range
blocksize = .2; %length of stimulus blocks in sec
nblock = round((blocksize/dtStim)); %number of blocks in simulation
stimrnge = [-25 25]; % range of stimulus y values
stim1 = reshape(ones(nblock,1)*rand(1,nt/nblock),[],1); % generate random stimulus in 1D
stim1 = stim1*diff(stimrnge)+stimrnge(1); % shift to have desired range
stim2 = reshape(ones(nblock,1)*rand(1,nt/nblock),[],1); % generate random stimulus in 1D
stim2 = stim2*diff(stimrnge)+stimrnge(1); % shift to have desired range
stim = [stim1-.5,-stim2+.5]'; % join both vectors into a 2-d stimulus

% plot stimulus
subplot(541);
h = plot(time, time*0, 'k--', time, stim);
set(h(2), 'linewidth', 2); box off
set(gca, 'FontSize',15, 'linewidth',1)
set(gca, 'xlim', [0 NT])
ylabel('stimulus')

%% Set up the population of GLM neurons
params.N = 400; % Number of neurons
params.taud = 1; % rate decay
params.mu = 0; % quadratic cost on spiking
params.tdel = 0; % time delay in bins
params.wsig = 2; % noise in weights 
% for local framework
params.fmax = 1; % max firing rate
params.alpha = 80; % alpha, controls precision
params.fmin = 0; % baseline firing rate 
% for population framework
params.beta = 50; % scaling of pspike

%% run simulations for both local and population frameworks
% 1-D integrator
params.A = 0; % integrator dynamics
params.wmean = 0.2*ones(1,size(params.A,1)); % mean weight value
[o, xh, x] = local_framework(stim(1,:), dtStim, params); %local framework simulation
[o1, xh1, x1] = population_framework(stim(1,:), dtStim, params); %population framework simulation

% 2-D harmonoic oscillator
params.A = inv([-.01 .1; -.1 -.01]); %oscillator dynamic
params.wmean = 0.2*ones(1,size(params.A,1)); %mean weight value
[o2, xh2, x2] = local_framework(stim, dtStim, params);% local framework simulation
[o3, xh3, x3] = population_framework(stim, dtStim, params);%population framework simulation

%% plot all four simulations - raster plots and read-outs
figure(1)
% raster plot for local framework 1-D integrator
subplot(5,4,[5 9]);
N = params.N; % number of neurons
[iiinh,jjinh] = find(o(1:N/2,:)); % spikes of inhibitory
[iiexc,jjexc] = find(o(N/2+1:N,:)); % spikes of excitatory
h = plot(jjinh*dtStim,iiinh, '.', jjexc*dtStim,iiexc+N/2,'.');
set(gca, 'ylim',[0 N])
box off
set(gca, 'FontSize',15, 'linewidth',1,'xlim', [0 NT])
ylabel('neuron')
title('local 1-D')

% read-out for local framework 1-D integrator 
subplot(5,4,[13 17]);
h1 = plot(time, x, 'k-');
hold on
h2 = plot(time, xh, '-',time, time*0, 'k--');
set(h1(1), 'linewidth', 1);
box off
set(gca, 'FontSize',15, 'linewidth',1, 'xlim', [0 NT], 'ylim',[min(xh(1,:))*1.2, max(xh(1,:))*1.2])
hold off
ylabel('read-out')
xlabel('time')

% raster plot for population framework 1-D integrator
subplot(5,4,[6 10]);
[iiinh,jjinh] = find(o1(1:N/2,:));
[iiexc,jjexc] = find(o1(N/2+1:N,:)); 
h = plot(jjinh*dtStim,iiinh, '.', jjexc*dtStim,iiexc+N/2,'.');
set(gca, 'ylim',[0 N])
box off
set(gca, 'FontSize',15, 'linewidth',1,'xlim', [0 NT])
title('population 1-D')

% read-out for population framework 1-D integrator
subplot(5,4,[14 18]);
h1 = plot(time, x1, 'k-');
hold on
h2 = plot(time, xh1, '-',time, time*0, 'k--');
set(h1(1), 'linewidth', 1);
box off
set(gca, 'FontSize',15, 'linewidth',1, 'xlim', [0 NT], 'ylim', [min(xh1(1,:))*1.2, max(xh1(1,:))*1.2])
hold off

% raster plot for local framework 2-D oscillator
subplot(5,4,[7 11]);
N = params.N;
[iiinh,jjinh] = find(o2(1:N/2,:)); % spikes of inhibitory
[iiexc,jjexc] = find(o2(N/2+1:N,:)); % spikes of excitatory
h = plot(jjinh*dtStim,iiinh, '.', jjexc*dtStim,iiexc+N/2,'.');
set(gca, 'ylim',[0 N])
box off
set(gca, 'FontSize',15, 'linewidth',1,'xlim', [0 NT])
title('local 2-D')

% read-out for local framework 2-D oscillator
subplot(5,4,[15 19]);
h1 = plot(time, x2, 'k-');
hold on
h2 = plot(time, xh2, '-',time, time*0, 'k--');
set(h1(1), 'linewidth', 1);
set(h1(2), 'linewidth', 1,'Color', 'k')
box off
set(gca, 'FontSize',15, 'linewidth',1, 'xlim', [0 NT], 'ylim',[min(xh2(1,:))*1.2, max(xh2(1,:))*1.2])
hold off

% raster plot for population framework 2-D oscillator
subplot(5,4,[8 12]);
[iiinh,jjinh] = find(o3(1:N/2,:));
[iiexc,jjexc] = find(o3(N/2+1:N,:)); 
h = plot(jjinh*dtStim,iiinh, '.', jjexc*dtStim,iiexc+N/2,'.');
set(gca, 'ylim',[0 N])
box off
set(gca, 'FontSize',15, 'linewidth',1,'xlim', [0 NT])
title('population 2-D')

% read-out for population framework 2-D oscillator
subplot(5,4,[16 20]);
h1 = plot(time, x3, 'k-');
hold on
h2 = plot(time, xh3, '-',time, time*0, 'k--');
set(h1(1), 'linewidth', 1);
set(h1(2), 'linewidth', 1,'Color', 'k')
box off
set(gca, 'FontSize',15, 'linewidth',1, 'xlim', [0 NT], 'ylim',[min(xh3(1,:))*1.2, max(xh3(1,:))*1.2])
hold off

%% Compute errors and spikes fired
err(1) = 1 - sum((x-xh).^2)/sum((xh - mean(xh).^2)); % error in local framework 1-D integrator
err(2) = 1 - sum((x1-xh1).^2)/sum((xh1 - mean(xh1)).^2);% error in population framework 1-D integrator
err(3) = 1 - sum((x2-xh2).^2)/sum((xh2 - mean(xh2)).^2);% error in local framework 2-D oscillator
err(4) = 1 - sum((x3-xh3).^2)/sum((xh3 - mean(xh3)).^2);% error in population framework 2-D oscillator

spcts = [sum(o(:)), sum(o1(:)),sum(o2(:)), sum(o3(:))];
fprintf('Errors:\nlocal 1D=%.3e, \npopulation 1D=%.3e\nlocal 2D=%.3e\npopulation 2D=%.3e', err);
fprintf('\nTotal spikes:\nlocal 1D =%d,\npopulation 1D =%d\nlocal 2D =%d\npopulation 2D=%d\n', spcts);

%% plot errors and spikes fired
figure(2)
subplot(121)
bar(err)
set(gca, 'linewidth', 1, 'FontSize', 15)
box off
title('R^2 value')
xticklabels({'local 1D','population 1D','local 2D','population 2D'})

subplot(122)
bar(spcts)
title('spikes fired')
set(gca, 'linewidth', 1, 'FontSize', 15)
xticklabels({'local 1D','population 1D','local 2D','population 2D'})
box off