function start

%%Parameters

set.lambda  = 0.1;       % influx rate of susceptible fly
set.gamma   = 0.2;       % rate of susceptible fly to becomes carrier
set.p       = 0.5;       % probability of fly that are female
set.psiE    = 0.4;       % average maturation rate from egg to larvae
set.psiL    = 0.6;       % average maturation rate from larvae to pupae
set.psiU    = 0.7;       % average maturation rate from pupae to adult fly
set.betaEf  = 0.1;       % natural death rate of egg fly
set.betaLf  = 0.1;       % natural death rate of larvae fly
set.betaUf  = 0.1;       % natural death rate of pupae fly
set.betaSf  = 0.1;       % natural death rate of susceptible fly
set.betaCf  = 0.3;       % oviposition rate of female fly
set.delta   = 0.3;       % oviposition rate of female fly
set.D1sf    = 0.001;     % diffusion coefficient among susceptible fly
set.D1cf    = 0.001;     % diffusion coefficient among carrier fly
set.xi      = 0.7;       % rate of carrier fly lay pathogen on water or food
set.alpha   = 0.03;       % efficiency of elimination of breeding site
set.phi = 0.1;           % efficiency of sanitation
set.vartheta = 0.1;      % efficiency of installation UV light trap

set.rho = 0.1;
set.tau = 0.1;           % efficiency of water purification
set.theta   = 0.1;       % influx rate of susceptible human
%set.mu      = 0.4;       % rate of susceptible human becomes infected
set.epsilon = 0.0008;       % rate of infected human becomes recovered
set.sigma   = 0.001;       % rate of recovered human becomes susceptible again
set.betaSH  = 0.0008;       % natural death rate of susceptible human
set.betaIH  = 0.0008;       % natural death rate of infected human
set.betaRH  = 0.0008;      % natural death rate of recovered human
set.D2sh    = 0.1;       % diffusion coefficient among susceptible human
set.D2ih    = 0.3;       % diffusion coefficient among infected human
set.D2rh    = 0.1;       % diffusion coefficient among recovered human
set.eta     = 0.55;           % rate of contaminated water or food to be consumed by susceptible human 
set.omegaIH = 0.3;       % disease-induced death rate of infected human

%%Initial conditions
initial.EfB = 100;        % set the initial value of 'EfB'
initial.LfB = 100;          % set the initial value of 'LfB'
initial.UfB = 10;          % set the initial value of 'UfB'
initial.Sf = 10;        % set the initial value of 'Sf'
initial.Cf = 10;         % set the initial value of 'Cf'

initial.ShB = 1000;       % set the initial value of 'Sh'
initial.IhB = 1;          % set the initial value of 'Ih'
initial.RhB = 0;          % set the initial value of 'Rh'

%%Set simulation time
end_time = 50;

%Definition of the ODE system

%%Function to calculate derivatives of the ELUSC-SIR model
function deriv = ode_system (t, x, set)

% Input:
%       t:     Time
%       x:     Vector of the current values of all variables in the same
%              order as defined the inital values: ELUSC-SIR
%       set: Used to pass parameter values.
% Output:
%       deriv: Column vector of derivatives, must be the same order as the
%              input vector x.

EfB = x(1);
LfB = x(2);
UfB = x(3);
Sf = x(4);
Cf = x(5);

ShB = x(6);
IhB = x(7); 
RhB = x(8);

%elimination of breeding site
dEfB = set.p * set.delta * Sf + set.p * set.delta * Cf - set.betaEf * EfB - set.psiE * EfB - set.alpha * EfB;
dLfB = set.psiE * EfB - set.betaLf * LfB - set.psiL * LfB - set.alpha * LfB;
dUfB = set.psiL * LfB - set.betaUf * UfB - set.psiU * UfB - set.alpha * UfB;
dSf = set.lambda * Sf + set.psiU *UfB - set.betaSf * Sf - set.gamma * Sf + set.D1sf * Sf;
dCf = set.gamma * Sf - set.betaCf * Cf + set.D1cf * Cf;

dShB = set.theta * ShB - set.betaSH * ShB - (set.xi * set.eta * ShB) + set.sigma * RhB + set.D2sh * ShB;
dIhB = set.xi * set.eta * ShB - set.betaIH * IhB - set.omegaIH * IhB - set.epsilon * IhB + set.D2ih * IhB; 
dRhB = set.epsilon * IhB - set.betaRH * RhB - set.sigma * RhB + set.D2rh * RhB; 

deriv = [dEfB; dLfB; dUfB; dSf; dCf; dShB; dIhB; dRhB];
end

% Solve the ODE system and plot the results

%N = initial.Ef + initial.Lf + initial.Uf + initial.Sf + initial.Cf + initial.Ef + initial.Lf + initial.Uf;
%R_0 = set.xi * set.eta / set.epsilon

% Extract initial values from the 'initial' structure and collect them
% in a column vector for use in 'ode45'.
initial_values = [];
variable_names = fieldnames(initial);
for i=1:length(variable_names) 
    initial_values = [initial_values; initial.(variable_names{i})];
end

% integrate the ODE system
[t, y] = ode45(@(t, x) ode_system(t, x, set), ...
               [0 end_time], ...
               initial_values, ...
               []);

% plot the results
plot(t, y);
xlabel('time');
ylabel('number of individuals');
title('Fly and Human system');
legend('Egg controlled','Larvae controlled','Pupae controlled','Susceptible fly',...
    'Carrier fly','Susceptible human controlled','Infected controlled','Recovered controlled');
hold on
end
