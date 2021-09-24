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

% Specify parameter set
model_params_yaml_filepath = fullfile('config','eg4_xflr.yaml');

% Construct model
% % kite = LonKiteDynamics(stateRep, model_params_yaml_filepath);
% % kite = LonKiteDynamicsWithPos(stateRep, model_params_yaml_filepath);
kite = KiteDynamics(model_params_yaml_filepath);

x0 = kite.defaultState; % This random state is feasible for all representations
u0 = kite.defaultControl;

% Test dynamics function
params = [1; 0; 0];
xdot = kite.dynamics(x0, u0, params)

%% Simulate model
sim_dt = 0.005; % 200 Hz
simulator = Simulator(sim_dt, @kite.dynamics);

% Test one step simulation
xp = simulator.simulate_step(x0, u0, params)
%xp = simulator.simulate_step(x0, u, @model.eval_dynamics)      % Override dynamics 
%xp = simulator.simulate_step(x0, u, @model.eval_dynamics, 0.1) % Override dynamics and dt

% Test simulation until end time
tf = 2.0;
T = 0:sim_dt:tf
U = repmat(u0, 1, length(T)-1)

X = simulator.simulate(x0, U, params);
%X = simulator.simulate(x0, U, @model.eval_dynamics)         % Override dynamics 
%X = simulator.simulate(x0, U, @model.eval_dynamics, dt)     % Override dynamics and dt

%%
% Plot trajectory (with/without plotting control input)
% simulator.plot_state_trajectory(T, X, kite.sys);
% simulator.plot_control_state_trajectory(T, U, X, kite.sys);
ph = kite.plot_pose_trajectory(X);
ph = kite.plot_pose_trajectory(X+2, '2', ph);

% Store all simulation data in common timetable
simulator.plot_state_trajectory(T, X, kite.sys);
simulator.plot_control_state_trajectory(T, U, X, kite.sys)
