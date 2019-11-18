%% Figure 5
%
% =======
% Characterizes the performance of the local framework. The first part of
% the code allows you to vary the parameters alpha, fmax, and fmin to
% observe qualitative changes in spiking and read-out quality. The second
% half produces two plots showing how decoding error and sparsity of
% activation vary with alpha and f_max.
% =======
% 
% Dependencies: local_framework.m, colorcet.m
%
%% Set up stimulus
NT = 1; % simulation length
dtStim = .1/1000; % dt 
time = (dtStim:dtStim:NT); % time vector
nt = length(time); % number of time bins

% create stimulus
stim = zeros(1, nt);
stim(.25/dtStim:.45/dtStim) = 25;
stim(.65/dtStim:.85/dtStim) = -50;

% plot stimulus
figure(1)
subplot(7,1,[1]);
h = plot(time, time*0, 'k--', time, stim);
set(h(2), 'linewidth', 2);
set(gca,'xticklabel', []);
set(gca, 'FontSize',15, 'linewidth',1)
box off;
ylabel('stimulus')

%% Set up the population of GLM neurons
% vary alpha, fmin or fmax to explore how these parameters affect spiking
% behavior and read-out quality

params.N = 400; % Number of neurons
params.wmean = 0.1; % mean weight value
params.wsig = .01; % variance in weight value
params.fmax = 100;% maximal firing rate 
params.alpha = 1000;% alpha, controls precision
params.fmin = 1; % baseline firing rate
params.taud = 10; % rate decay
params.A = 1/Inf; % dynamics of x 
params.mu = 0; % quadratic cost on spiking
params.tdel = 0; % time delay, in bins

%% Simulate GLM population with parameters set above
[o, xh, x] = local_framework(stim, dtStim, params);

%% Plot raster plot and read-out

%raster plot
subplot(7,1,[2 3 4])
N = params.N; % number of neurons
[iiinh,jjinh] = find(o(1:N/2,:)); % spikes of inhibitory
[iiexc,jjexc] = find(o(N/2+1:N,:)); % spikes of excitatory
plot(jjinh*dtStim,iiinh, '.', jjexc*dtStim,iiexc+N/2,'.');
set(gca, 'ylim',[0 N])
box off
set(gca, 'FontSize',15, 'linewidth',1)
ylabel('neuron')

% target and read-out
subplot(7,1,[5 6 7])
h = plot(time, x, 'k-', time, xh, '-', time, time*0, 'k--');
legend('target', 'read-out')
legend boxoff
box off
set(gca, 'FontSize',15, 'linewidth',1)
ylim([-6 6])
set(gcf, 'renderer', 'Painters')
xlabel('time')

%% Generate error plots
% cycles through a range of values for alpha and f_max to generate two
% plots, one showing how read-out error varies with these values, and
% another showing how population activity varies.

astep = 500; % step size in alpha values
alphavals = 10:astep:10^4; 
fstep = 500; % step size in f_max values
fvals = 10:fstep:10^4; 

net_err = zeros(length(alphavals), length(fvals)); % net read-out error 
perc_active = zeros(length(alphavals), length(fvals)); % percent of population active

for i = 1:length(alphavals) %cycle through alpha values
    params.alpha = alphavals(i); %set simulation parameter to this alpha value
for j = 1:length(fvals) %cycle through fmax values
   params.fmax = fvals(j); %set simulation parameter to this fmax value
   [o, xh, x] = local_framework(stim, dtStim, params); %run simulation with these parameters
   net_err(i,j) = 1 - sum((xh-x).^2)/sum((xh - mean(xh)).^2); %get decoding error with these params
   perc_active(i,j) = mean(sum(o))/N; %get percent active with these params 
end
end

%% plot decoding error
figure(2)
title('decoding error')
pcolor(alphavals, fvals,net_err);
shading interp;
hold on
scatter(1000,100,'k*')

set(gca,'XScale','log', 'YScale', 'log','FontSize', 15)
colorcet('L12', 'reverse', 1); %perceptually uniform colormap
colorbar;
set(gcf,'renderer','Painters')

%% plot % of population active 
figure(3)
title('% of population active')
pcolor(alphavals, fvals, perc_active)
shading interp;
hold on
scatter(1000,100,'k*')
set(gca,'XScale','log', 'YScale', 'log','FontSize', 15)
set(gca, 'clim', [0 max(max(perc_active))]);
colorcet('L12')
set(gcf,'renderer','Painters')
colorbar;