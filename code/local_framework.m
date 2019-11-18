function [ss,xh,xx] = local_framework(stim, dt, params)
% [ss,xh,xx] = local_framework(stim,dt, params)
%
% Function simulates network of N GLM neurons tracking a linear
% transformation of a D dimensional target, x, which is a function of the
% time varying stimulus, c. The function outputs the values of x, the
% network read-out, xh, and the spike trains, s. Spiking is
% conditioned on the individual neuron's internal copy of the estimation
% error.
%
% Inputs:
% =======
% stim [D x T] - stimulus (# inputs x # time points)
% dt - time step size
% params - struct of simulation parameters with fields:
%  .N - number of neurons
%  .alpha - constant modulating 'hardness' of threshold
%  .fmax - maximal firing rate per neuron (in spikes/unit time)
%  .fmin - background firing rate (spikes/unit time)
%  .taud - decay rate of filtered spike train
%  .A - dynamics of x
%  .mu - quadratic cost term
%  .wmean - mean weight value
%  .wsig - stdev of Gaussian noise added to weights
%  .tdel - time delay (in bins)
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

% Params governing GLM nonlinearity
alpha = params.alpha; % gain for nonlinearity
fmax = params.fmax; % maximal firing rate for nonlinearity
fmin = params.fmin; % background spiking

% dynamics parameters 
lambda_d = 1/params.taud; % decay in filtered spike train
A = params.A; % dynamics matrix of x

% cost params
mu = params.mu; %quadratic cost term
rlambda = params.mu*lambda_d^2; %cost on spiking

%% Generate decoding weights 
ww = params.wmean'*[-ones(1,ceil(N/2)), ones(1,floor(N/2))] ... % mean of weights
    + randn(nd,N)*params.wsig; % noise

T = (mu*lambda_d^2 +sum(ww'.^2,2))/2; %Threshold

%% Initialize variables
ss = zeros(N,nt); % spikes
rr = zeros(N,nt); % filtered spike train
xx = zeros(nd,nt); % true x
xh = zeros(nd,nt); % network output (estimate of x)
zz = zeros(nd,nt); % internal estimate of x 
vv = zeros(N, nt); % voltage
pspike = zeros(N,nt); % probability of spiking
pen = zeros(N,1); % penalty on spking 

% exponential integrators
Amult = expm(dt*A); % x integration 
Amult_delay = expm(tdel*dt*A); % x integration with delay 
Rmult = exp(-lambda_d*dt); % r integration
Rmult_delay = exp(-lambda_d*dt*tdel); % r integration with delay

%% Run simulation

% Loop over time bins
for tt = 2:nt
    
    % === Update true target variable ========    
    xx(:,tt) = Amult*xx(:,tt-1)+dt*stim(:,tt);
    
    % === Update filtered spike trains ======
    rr(:,tt) = rr(:,tt) + Rmult*rr(:,tt-1);
    xh(:,tt) = ww*rr(:,tt); % current network estimate (pre-spike)
        
    % === Update z(t) (network estimate of x(t)) =====
    zz(:,tt) = Amult*zz(:,tt-1)+dt*stim(:,tt);

    % === Update penalty =====
    dpen = - rlambda*(ss(:,tt-1)+Rmult_delay*rr(:,tt-1)); % deriv of penalty
    pen = pen+dt*dpen;
    
    % === Update voltage ===== 
    vv(:,tt) = ww'*(Amult_delay*zz(:,tt)-Rmult_delay*ww*rr(:,tt))+pen;
    
    % === Compute instantaneous spike rate =====
    rt = alpha*(vv(:,tt) - T); 
    cond = fmax*1./(1+fmax*exp(-rt)) + fmin; %conditional intensity
    pspike(:,tt) = 1-exp(-cond*dt); %calculate probability of spiking

    % === Spiking =====
    iisp = find(rand(N,1)<pspike(:,tt)); %find neurons that spiked
    if ~isempty(iisp) 
        if tt < nt - tdel % avoid indexing error
        ss(iisp, tt+tdel) = 1; % insert spikes
        rr(iisp, tt+tdel) = rr(iisp, tt+tdel) + 1; % update 
        end
    end
        
    % === Update x estimate=====
    xh(:,tt) = ww*rr(:,tt);
end
