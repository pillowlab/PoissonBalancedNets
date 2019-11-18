%% Figure 4 
%
% =======
% Compares the BSN LIF to the GLM single neuron. The neuron
% performing an integration of an OU process. By varying the parameter
% alpha, you can observe how the spike spread (and by extension, coding 
% accuracy) changes with respect to the hardness of the GLM threshold. 
% =======
% 
% Dependencies: plotRaster
%
%% Set up simulation 
NT = 15; %simulation length
dt = .1; %dt 
time = dt:dt:NT; %time vector
T = 0.2; % threshold
nt = length(time); % number of time bins

c = zeros(1,nt); %stimulus
lambda_c = 1/20; % decay of stimulus

% --- Make OU stimulus -----
for i = 2:nt
    dc = -lambda_c*c(i-1) + 0.1*randn(); %derivative 
    c(i) = c(i-1) + dc; %update c
end

%plot stimulus
clf
subplot(511);
plot(c);
box off
title('stimulus')

%% Simulate LIF neuron
lif_v = zeros(1,nt); % voltages
lif_sp = zeros(1,nt); % spike times

lambda_v = 0.05; % voltage decay

for i = 2:nt

    dv = -lambda_v*lif_v(i-1) + c(i); % derivative of voltage
    lif_v(i) = lif_v(i-1) + dv*dt; %integration of voltage
    
    if lif_v(i) > T % check if neuron is above threshold
        lif_sp(i) = 1; % record spike 
        lif_v(i) = -2*T; % reset to -2*Threshold
    end
end

tsp = dt*find(lif_sp); % spike times

% Plot LIF voltage trace and spikes
subplot(512);
plot(time, lif_v, [0 NT], [T T], 'k--')
hold on
plotRaster(tsp)
box off
title('LIF')
    
%% Simulate GLM 
alpha_val = 100; % VARY THIS to change the hardness of the spiking threshold
fmax = 50; % upper limit on spiking
fmin = 0; % lower limit on spiking
reps = 5; % repetitions of GLM simulation

sptimes = cell(1,reps); % cell array with spike times for each rep (for raster plot)
vvals = cell(1,reps); % cell array with voltage traces for each repetition

for j = 1:reps
    % initialize for each rep
    glm_v = zeros(1,nt); % Voltage
    glm_sp = zeros(1,nt); % Spikes
    
    for i = 2:nt
        dv = -lambda_v*glm_v(i-1) + c(i); % derivative of voltage
        glm_v(i) = glm_v(i-1) + dv*dt; % integrate voltage
        
        rt = alpha_val*(glm_v(i) - T); %instantaneous sp rate 
        cond = fmax*1./(1+fmax*exp(-rt)) + fmin; %conditional intensity
        Ps = 1-exp(-cond*dt); %probability of spiking

        if rand() < Ps % threshold
            glm_sp(i) = 1; %insert spike
            glm_v(i) = -2*T; %voltage reset
        end
        
    end
    
    % Store spikes
    sptimes{j} = find(glm_sp)*dt;
    vvals{j} = glm_v;
end

% plot GLM voltage trace and spikes for all reps
subplot(5,1,[ 3 4 5]);
for i = 1:reps
plot(time, vvals{i}, [0 NT], [T T], 'k--')
hold on
plotRaster(sptimes{i}, 'k', [], [], i)
end
box off
set(gca, 'ylim', [min(lif_v) reps + 0.5])
title('GLM')
