function [ss,xh,xx] = population_framework(stim, dt, params)
% [ss, xh, xx] = population_framework(stim, dt, params)
%
% Function simulates network of N GLM neurons tracking a linear
% transformation of an D-dimensional target, x, which is a function of the
% time varying stimulus, c. The function outputs the values of x, the
% network read-out, xh, and the spike trains, s. This approach exploits the
% pseudoinverse of weights and approaches spiking from a population level. 
%
% Inputs:
% =======
% stim [D x T] - stimulus (# inputs x # time points)
% dt - time step size
% params - struct of simulation parameters with fields:
%  .N - number of neurons
%  .taud - decay rate of filtered spike train
%  .A - dynamics of x
%  .wmean - mean weight value
%  .wsig - stdev of Gaussian noise added to weights
%  .beta - timescale of error correction in bins
%  .tdel - time delay in bins
%
% Outputs:
% ========
% ss [N x T] - spike trains for all N neurons, length time
% xh [D x T] - network output (reconstruction of x)
% xx [D x T] - true dynamic variable x

%% Extract parameters
[nd,nt] = size(stim); % number of dimensions and time bins for stimulus
N = params.N; % number of neurons
tdel = params.tdel; % time delay in bins
beta = params.beta; % error timescale in bins

% dynamics parameters 
lambda_d = 1/params.taud; % decay in filtered spikes
A = params.A; % dynamics matrix of x(t)

%% Generate decoding weights 
ww = params.wmean'*ones(1,N/2) ... % mean of weights
    + randn(nd, N/2)*params.wsig*mean(params.wmean); % noise

wpinv = pinv(ww);  % pseudoinverse of weights

%% Initialize variables

ss = zeros(N,nt,1); % spikes, includes the mirror neurons
rr = zeros(N/2,nt); % filtered spike train
vv = zeros(N/2,nt); % voltages
zz = zeros(nd,nt); % network estimate of x
xh = zeros(nd, nt); % estimate
xx = zeros(nd,nt); % true x
pspike = zeros(N/2,nt); % probability of spiking

% exponential integrators
Amult = expm(dt*A); % x integration 
Amult_delay = expm(tdel*dt*A); % x integration with delay
Rmult = exp(-lambda_d*dt); % r integration
Rmult_delay = exp(-lambda_d*dt*tdel); % r integration with delay

%% Run simulation

for tt = 2:nt
    
    % === Update true target variable ========    
    xx(:,tt) = Amult*xx(:,tt-1)+dt*stim(:,tt);
    
    % === Update filtered spike trains ======
    rr(:,tt) = rr(:,tt) + Rmult*rr(:,tt-1); %integrate rate
    xh(:,tt) =  ww*rr(:,tt); % current network estimate (if there are no spikes)
    
    % === Update z(t) (network estimate of x(t)) =====
    zz(:,tt) = Amult*zz(:,tt-1)+dt*stim(:,tt);
    
    % === Voltage update ===== 
    vv(:,tt) = wpinv*(Amult_delay*zz(:,tt)-Rmult_delay*ww*rr(:,tt));    
    
    % === Compute spike probabilities ====
    pspike(:,tt) = (1/beta)*vv(:,tt); %probability of spiking
    ppos = max(pspike(:,tt),0); % p(spike) for "pos" mirror neurons
    pneg = -min(pspike(:,tt),0); % p(spike) for "neg" mirror neurons
    
    % === Sample spikes ====
    sspos = rand(N/2,1)<ppos; % positive mirror neuron spikes
    ssneg = rand(N/2,1)<pneg; % negative mirror neuron spikes
    if tt < nt - tdel % avoid index out of bounds
    ss(:,tt + tdel) = [ssneg;sspos]; % add sps to spike train in future
    netsps = sspos-ssneg; %net # of spikes
    rr(:, tt + tdel) = rr(:, tt + tdel) + netsps;% update filtered spike train tdel in future
    end

     % === Update x estimate using filtered spike trains ==== 
    xh(:,tt) = ww*rr(:,tt);

end
    

