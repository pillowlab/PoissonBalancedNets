function [ss, xh, xx] = bsn(stim, params)
% [ss,xh, xx] = bsn(stim,dt, params)
%
% Function simulates balanced spiking network of N LIF neurons integrating
% a D dimensional target, x, which is a function of the
% time varying stimulus, c. The function outputs the values of x, the
% network read-out xh, and the spike trains, s. This code is based on the 
% model presented by Boerlin et al. in 'Predictive coding of dynamical
% variables in balanced spiking networks' (2013). 
% 
% Inputs:
% =======
% stim [D x T] - stimulus (# inputs x # time points)
% dt - time step size
% params - struct of simulation parameters with fields:
%  .N - number of neurons
%  .mu - quadratic cost term
%  .nu - linear cost term
%  .wmean - mean weight value
%  .wsig - stdev of Gaussian noise added to weights
%  .taud - spike train time decay constant
%  .tau_v - voltage time decay constant
%  .sig_v - voltage noise
%  .sp_cond - spiking condition 
%  
% Outputs:
% ========
% ss [N x T] - spike trains for all N neurons, length time
% xh [D x T] - network output
% xx [D x T] - target
%% Get parameters
dt= params.dt; %time bin size
[nd, nt] = size(stim); % number of dimensions and time bins
N = params.N; % number of neurrons

lambda_d = 1/params.taud; % decay in filtered spike train
lambda_v = 1/params.tau_v; %decay in voltage

nu = params.nu; % linear cost term
mu = params.mu; % quadratic cost term

sig_v = params.sigv; % voltage noise

spiking_condition = params.sp_cond;  % set spike condition (random selection or all)

%% Generate decoding weights
ww = params.wmean'*[-ones(1,ceil(N/2)), ones(1,floor(N/2))] ... %mean
    + randn(nd, N)*params.wsig;
I = speye(N); %identity matrix of size N
Thresh =  (mu*lambda_d^2 + nu*lambda_d + sum(ww'.^2,2))/2; %threshold
%%
% initialize variables
Vt = zeros(N,nt); % voltage
ss = zeros(N,nt); % spikes
rt = zeros(N,nt); % spike rate (filtered spikes)
xx = zeros(nd, nt); % target

% initialize target
xx(:,1) = stim(:,1)*dt;
%% Run simulation
for tt = 2:nt

    % === Update true target variable ========
    xx(:,tt)  = xx(:, tt-1) + stim(:,tt)*dt;
    
    % === Update filtered spike trains ======
    rt(:,tt) = I*(1-lambda_d*dt)*rt(:,tt-1);
    
    % === Update voltage ======
    % Leak, input, slow dynamics and noise
    dV = -lambda_v*dt*Vt(:,tt-1) + ww'*dt*stim(:,tt) + ...
        (ww'*ww)*lambda_d*dt*rt(:,tt-1) + sig_v*sqrt(dt)*randn(N,1); %derivative
    Vt(:,tt) = Vt(:,tt-1) + dV;

    % === Spiking =====
    iisp = find(Vt(:,tt)>Thresh(:)); % find neurons over threshold 
    
    % === Spiking conditions =====
    if ~isempty(iisp)
        if spiking_condition == 1  %Randomly select from neurons over T
            spikeidx = randsample(iisp,1);
            ss(spikeidx,tt) = 1; % insert spike
            rt(spikeidx,tt) = rt(spikeidx,tt)+1;  % update rates
            Vt(:,tt) = Vt(:,tt) - (ww'*ww + mu*lambda_d^2*I)*ss(:,tt); % membrane reset 
        else %let all fire
            spikeidx = iisp;
            ss(spikeidx,tt) = 1; % insert spike
            rt(spikeidx,tt) = rt(spikeidx,tt)+1;  % update rates
            Vt(:,tt) = Vt(:,tt) - (ww'*ww + mu*lambda_d^2*I)*ss(:,tt); % membrane reset
        end
    end

end
% get read-out
xh = ww*rt;

end