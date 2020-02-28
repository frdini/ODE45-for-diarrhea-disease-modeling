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
set.betaCf  = 0.1;       % natural death rate of carrier fly
set.delta   = 0.3;       % oviposition rate of female fly
set.D1sf    = 0.001;     % diffusion coefficient among susceptible fly
set.D1cf    = 0.001;     % diffusion coefficient among carrier fly
set.xi      = 0.7;       % rate of carrier fly lay pathogen on water or food
set.alpha   = 0.03;       % efficiency of elimination of breeding site
set.phi = 0.03;           % efficiency of sanitation
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
set.eta     = 0.6;           % rate of contaminated water or food to be consumed by susceptible human 
set.omegaIH = 0.3;       % disease-induced death rate of infected human

%Initial conditions
initial.Ef = 100;        % set the initial value of 'Ef'
initial.Lf = 100;          % set the initial value of 'Lf'
initial.Uf = 10;          % set the initial value of 'Uf'
initial.SfS = 10;        % set the initial value of 'EfB'
initial.CfS = 10;          % set the initial value of 'LfB'
initial.ShS = 1000;       % set the initial value of 'Sh'
initial.IhS = 1;          % set the initial value of 'Ih'
initial.RhS = 0;          % set the initial value of 'Rh'


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

Ef = x(1);
Lf = x(2);
Uf = x(3);
SfS = x(4);
CfS = x(5);
ShS = x(6);
IhS = x(7);
RhS = x(8);


dEf = set.p * set.delta * SfS + set.p * set.delta * CfS - set.betaEf * Ef - set.psiE * Ef;
dLf = set.psiE * Ef - set.betaLf * Lf - set.psiL * Lf;
dUf = set.psiL * Lf - set.betaUf * Uf - set.psiU * Uf;
%sanitation
dSfS = set.lambda * SfS + set.psiU *Uf - set.betaSf * SfS - set.gamma * SfS + set.D1sf * SfS - set.phi * SfS;
dCfS = set.gamma * SfS - set.betaCf * CfS + set.D1cf * CfS - set.phi * CfS;

%sanitation
dShS = set.theta * ShS - set.betaSH * ShS - (set.xi * set.eta * ShS) + set.sigma * RhS + set.D2sh * ShS;
dIhS = set.xi * set.eta * ShS - set.betaIH * IhS - set.omegaIH * IhS - set.epsilon * IhS + set.D2ih * IhS; 
dRhS = set.epsilon * IhS - set.betaRH * RhS - set.sigma * RhS + set.D2rh * RhS; 

deriv = [dEf; dLf; dUf; dSfS; dCfS; dShS; dIhS; dRhS];
end

% Solve the ODE system and plot the results

%N = initial.Sf + initial.Cf + initial.SfS + initial.CfS;
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
legend('Egg','Larvae','Pupae','Susceptible Flies controlled','Carrier Flies controlled',...
    'Susceptible Human controlled','Infected controlled','Recovered controlled')
hold on
end
