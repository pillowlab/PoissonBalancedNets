%% Evaluate Euler approximation
%
% =======
% Calculates the theoretical minimum error that arises from integrating x
% using exponential Euler integration with a given time delay.
% =======
%

function [xx,xxt] = test_approx(stim, dt, params)
A = params.A; %dynamics of x
[nd,nt] = size(stim); % number of dimensions of stimulus and length of simulation
tdel = params.tdel; %time delay

%Initialize vars
xx = zeros(nd,nt); % x with no time delay
xxt = zeros(nd,nt); % x with time delay
dxt = zeros(nd,nt); % derivative of x with delay

for tt = 2:nt-tdel
    % no time delay
    dx = A*xx(:,tt-1) +  stim(:,tt-1); %derivative 
    xx(:,tt) = xx(:,tt-1) + dx*dt; % intevrate
    
    %time delay
    dxt(:,tt+tdel) = expm(dt*tdel*A)*(A*xx(:,tt-1) + stim(:,tt-1)); % derivative at a future time
    xxt(:,tt) = xxt(:,tt-1) + dxt(:,tt)*dt; %integrate using derivative at current time
end 

end
