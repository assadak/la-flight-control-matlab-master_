clear all
close all
clc

addpath(fullfile('3rd_party', 'YAMLMatlab_0.4.3'))
addpath(genpath('src'))

%% Model initialization
% Choose state representation
%stateRep = 'UW';            % (vx, vz, wy, pitch)
%stateRep = 'Pitch';         % (Va, aoa, wy, pitch)
stateRep = 'Flightpath';    % (Va, aoa, wy, flightpath)

% Construct model(s)
% kite_xflr = LonKiteDynamicsWithPos(stateRep, fullfile('config','eg4_xflr.yaml'));
% kite_best = LonKiteDynamicsWithPos(stateRep, fullfile('config','eg4_xflr-Pvw-YR.yaml'));
kite_xflr = KiteDynamics(fullfile('config','eg4_xflr.yaml'));
kite_best = KiteDynamics(fullfile('config','eg4_xflr-Pvw-YR.yaml'));

x0 = kite_xflr.defaultState; % This random state is feasible for all representations
u0 = kite_xflr.defaultControl; % Neutral
params = [0;0;0];

%% Simulate model
sim_dt = 0.01;
simulator = Simulator(sim_dt);

% Test one step simulation
xp_xflr = simulator.simulate_step(x0, u0, params, @kite_xflr.dynamics);      % Override dynamics
xp_best = simulator.simulate_step(x0, u0, params, @kite_best.dynamics);

% Test simulation until end time
tf = 5.0;
T = 0:sim_dt:tf;

U = repmat(u0, 1, length(T)-1);     % Constant input

ident_dt = 0.16; % dt in 3*dt, 2*dt, 1*dt, 1*dt
ident2211 = kron([1 1 -1 -1 1 -1], [ones(1, ident_dt/sim_dt)]);
ident3211 = kron([1 1 1 -1 -1 1 -1], [ones(1, ident_dt/sim_dt)]);
identSignal = ident3211;
U(1, 1:size(identSignal,2)) = 0.3665 * 0.5 * identSignal;       % Max deflection = 0.3665 rad, choose 0.5 of max deflection as identification input

X_xflr = simulator.simulate(x0, U, params, @kite_xflr.dynamics);         % Override dynamics
X_best = simulator.simulate(x0, U, params, @kite_best.dynamics);

% Plot trajectories
simulator.plot_state_trajectory(T, {{X_xflr, 'xflr'}, {X_best, 'best'}}, kite_xflr.sys);
simulator.plot_control_state_trajectory(T, U, {{X_xflr, 'xflr'}, {X_best, 'best'}}, kite_xflr.sys);

ph = kite_xflr.plot_pose_trajectory(X_xflr, 'xflr');
ph = kite_xflr.plot_pose_trajectory(X_best, 'best', ph);

