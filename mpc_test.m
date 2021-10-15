clear all
close all
clc

addpath(fullfile('3rd_party', 'YAMLMatlab_0.4.3')) %% just adding path to needed functions
addpath(genpath('src'))

%% Model initialization
% Choose state representation
% % stateRep = 'UW';            % (vx, vz, wy, pitch)
% % stateRep = 'Pitch';         % (Va, aoa, wy, pitch)
stateRep = 'Flightpath';    % (Va, aoa, wy, flightpath)

% Construct control model
ctrlKite = LonKiteDynamics(stateRep, fullfile('config','eg4_xflr.yaml'));
% % ctrlKite = LonKiteDynamicsWithPos(stateRep, fullfile('config','eg4_xflr.yaml'));

% Construct simulation model
% % simKite = LonKiteDynamicsWithPos(stateRep, fullfile('config','eg4_xflr.yaml'));        % Most basic model
simKite = LonKiteDynamicsWithPos(stateRep, fullfile('config','eg4_xflr-Pvw-YR.yaml')); % Best identified

clear stateRep
%% Controller settings presets
ctrlset.normal.H          = 1.0; %horizon in time
ctrlset.normal.dt         = 0.05; % time steps
ctrlset.normal.N          = ceil( ctrlset.normal.H/ctrlset.normal.dt ); % how many steps 
ctrlset.normal.Nseg       = 2;
ctrlset.normal.Dpoly      = 5;
ctrlset.normal.comp_delay = 0.035;

ctrlset.nodelay = ctrlset; ctrlset.nodelay.comp_delay = 0;
ctrlset.shortH  = ctrlset; ctrlset.shortH.H = 0.2;
ctrlset.longH   = ctrlset; ctrlset.longH.H  = 4;

%% Controller initialization
% Choose controller settings preset here
ctrlset = ctrlset.normal;
% % ctrlset = ctrlset.nodelay;
% % ctrlset = ctrlset.shortH;
% % ctrlset = ctrlset.longH;

% Control objectives
ctrlobj = LonKiteOcp.ControlObjectives;
ctrlobj.track_angle      = 0; % (conflicts with tracking airspeed)
ctrlobj.track_airspeed   = 1;
ctrlobj.minimize_horizon = 1;
ctrlobj.track_height     = 0; % (only available with pos model)

% Construct controller(s)
opts.N    = ctrlset.N;
opts.Nseg = ctrlset.Nseg;
opts.D    = ctrlset.Dpoly;

controllers = {
    LonKiteOcp(ctrlKite, ctrlobj, ctrlset.H, opts, 'transcription','multiple_shooting');
    LonKiteOcp(ctrlKite, ctrlobj, ctrlset.H, opts, 'transcription','multiple_shooting', 'use_cos_grid',true, 'N',ceil(ctrlset.N/3));
    LonKiteOcp(ctrlKite, ctrlobj, ctrlset.H, opts, 'transcription','collocation', 'approx_control', true);
    };

% Set config and tuning parameters
ref.angle = deg2rad(-20);
ref.h     = 95; % Only with pos model
for iCtrl = 1:length(controllers)
    % Set reference
    controllers{iCtrl}.set_angle_ref( ref.angle );
    controllers{iCtrl}.set_Va_ref( 15 );
    controllers{iCtrl}.set_h_ref( ref.h );
    
    % Set tuning weights
    controllers{iCtrl}.set_mayer_multiplier( 10 );
    
    controllers{iCtrl}.set_W_angle_err( 10 );
    controllers{iCtrl}.set_W_Va_err( 1 );
    controllers{iCtrl}.set_W_h_err( 1 );
    
    controllers{iCtrl}.set_R_diag( [1, 0.1] );
end, clear i

% Set state and control bounds
ubx = [30, deg2rad(20), deg2rad(80), deg2rad(50)];
lbx = [9, deg2rad(-20), deg2rad(-80), deg2rad(-50)];
% % slack_bounds = [5, deg2rad(10), deg2rad(30), deg2rad(10)];
slack_bounds = [0, 0, 0, 0];

for iCtrl = 1:length(controllers)
    % Set state bounds with slacks. Slacks are set for nonzero entries in slack_bounds.
    % %     controllers{iCtrl} = controllers{iCtrl}.set_ubx( ubx, slack_bounds );
    % %     controllers{iCtrl} = controllers{iCtrl}.set_lbx( lbx, slack_bounds );
    
    controllers{iCtrl}.set_ubu( [ctrlKite.phyUBU(1), 0] );
    controllers{iCtrl}.set_lbu( ctrlKite.phyLBU );
end, clear i ubx lbx slack_bounds

% Set initial guess for state trajectory (longitudinal velocity needs to be nonzero, otherwise model produces nans)
x0 = ctrlKite.defaultState;
% % x0(4) = deg2rad(-60);
controllers{1}.set_X_guess( repmat(x0, 1, controllers{1}.N+1) );
controllers{2}.set_X_guess( repmat(x0, 1, controllers{2}.N+1) );
controllers{3}.set_X_guess( repmat(x0, 1, controllers{3}.N) );

if true
  %% Test single controller evalution (Run this section again for warm start)
    % Note:  - Resampled trajs are not plotted afterwards, so dt is chosen to be quasi-continuous.
    %        - For piecewise constant controls, duplicate last entry for better plotting.
    nRuns = 2;
    clear single
    for jRun = 1:nRuns
        for iCtrl = 1:length(controllers)
            controllers{iCtrl}.set_x0( x0 );
            t0 = tic;
            controllers{iCtrl} = controllers{iCtrl}.solve();
            single{iCtrl}.comp_time{jRun} = toc(t0);
        end
        dt = 0.0001; % dt for resampling trajectories
        % Multiple shooting
        [~, ~, Uopt, Topt] = controllers{1}.get_control_trajectory( dt );
        single{1}.U{jRun} = array2timetable([Uopt Uopt(:,end)]', 'RowTimes', seconds([Topt ctrlset.H])', 'VariableNames', ctrlKite.sys.InputName);
        [~, ~, Xopt, Topt] = controllers{1}.get_state_trajectory( dt );
        single{1}.X{jRun} = array2timetable(Xopt', 'RowTimes', seconds(Topt)', 'VariableNames', ctrlKite.sys.StateName);
        % Multiple shooting with cosine grid
        [U, T, Uopt, Topt] = controllers{2}.get_control_trajectory( dt );
        single{2}.U{jRun}    = array2timetable(U',    'RowTimes', seconds(T)',    'VariableNames', ctrlKite.sys.InputName);
        single{2}.Ucos{jRun} = array2timetable([Uopt Uopt(:,end)]', 'RowTimes', seconds([Topt ctrlset.H])', 'VariableNames', ctrlKite.sys.InputName);
        [X, T, Xopt, Topt] = controllers{2}.get_state_trajectory( dt );
        single{2}.X{jRun}    = array2timetable(X',    'RowTimes', seconds(T)',    'VariableNames', ctrlKite.sys.StateName);
        single{2}.Xcos{jRun} = array2timetable(Xopt', 'RowTimes', seconds(Topt)', 'VariableNames', ctrlKite.sys.StateName);
        % Collocation
        [U, T, Uopt, Topt] = controllers{3}.get_control_trajectory( dt );
        single{3}.U{jRun}    = array2timetable(U',    'RowTimes', seconds(T)',    'VariableNames', ctrlKite.sys.InputName);
        single{3}.Ucol{jRun} = array2timetable([Uopt Uopt(:,end)]', 'RowTimes', seconds([Topt ctrlset.H])', 'VariableNames', ctrlKite.sys.InputName);
        [X, T, Xopt, Topt] = controllers{3}.get_state_trajectory( dt );
        single{3}.X{jRun}    = array2timetable(X',    'RowTimes', seconds(T)',    'VariableNames', ctrlKite.sys.StateName);
        single{3}.Xcol{jRun} = array2timetable(Xopt', 'RowTimes', seconds(Topt)', 'VariableNames', ctrlKite.sys.StateName);
    end, clear jRun iCtrl t0 T X U Topt Xopt Uopt dt
    
    single_path = fullfile('cpp_test','single.yaml');
    if isfile(single_path)
        % If C++ test available, load trajectories
        fprintf(['Loading ' single_path '...'])
        single_yaml = ReadYaml(single_path);
        fprintf(' Done\n')
        TUopt = cell2mat(single_yaml.TU_opt)';
        Topt = TUopt(1,:);
        Uopt = TUopt(2:end,:);
        single_cpp.U = array2timetable([Uopt]', 'RowTimes', seconds([Topt])', 'VariableNames', ctrlKite.sys.InputName);
        
        TXopt = cell2mat(single_yaml.TX_opt)';
        Topt = TXopt(1,:);
        Xopt = TXopt(2:end,:);
        single_cpp.X = array2timetable(Xopt', 'RowTimes', seconds(Topt)', 'VariableNames', ctrlKite.sys.StateName);
        
        clear TUopt Topt Uopt U T
    end
    
    % Plot
    runToPlot = 2;
    %close all
    % Control trajectories ================================================
    figure('Name', 'Optimal control trajectory')
    subplot(1,1,1);
    for jRun = runToPlot:nRuns
        % Multiple shooting
        stairs(single{1}.U{jRun}.Time, single{1}.U{jRun}.dE, 'o-', 'MarkerSize', 6, 'LineWidth', 1, ...
            'DisplayName', ['mshoot it ' num2str(jRun) ' (' num2str(round(single{1}.comp_time{jRun} * 1000)) 'ms)']); hold on
        % Multiple shooting with cosine grid
        % %         stairs(single{2}.Ucos{jRun}.Time, single{2}.Ucos{jRun}.dE, 'o-', 'MarkerSize', 6, 'LineWidth', 1, ...
        % %             'DisplayName', ['ms\_cos it ' num2str(jRun) ' (' num2str(round(single{2}.comp_time{jRun} * 1000)) 'ms)']), hold on
        % Collocation
        % %         if controllers{3}.opts.approx_control
        % %             plot(single{3}.U{jRun}.Time, single{3}.U{jRun}.dE, 'LineWidth', 2, ...
        % %                 'DisplayName', ['colloc it ' num2str(jRun) ' (' num2str(round(single{3}.comp_time{jRun} * 1000)) 'ms)']), hold on
        % %         else
        % %             stairs(single{3}.Ucol{jRun}.Time, single{3}.Ucol{jRun}.dE, 'LineWidth', 2, ...
        % %                 'DisplayName', ['colloc it ' num2str(jRun) ' (' num2str(round(single{3}.comp_time{jRun} * 1000)) 'ms)']), hold on
        % %         end
        % %         plot(single{3}.Ucol{jRun}.Time, single{3}.Ucol{jRun}.dE, '^k', 'MarkerSize', 6, 'LineWidth', 1, 'DisplayName', 'colloc grid')
        
        if exist('single_cpp')
            % Multiple shooting
            stairs(single_cpp.U.Time, single_cpp.U.dE, 'o-', 'MarkerSize', 6, 'LineWidth', 1, ...
                'DisplayName', ['cpp mshoot it ' num2str(jRun)]); hold on
        end
        
    end, clear jRun
    grid on, legend('-DynamicLegend')
    
    % State trajectories ==================================================
    figure('Name', 'Optimal state trajectory')
    jRun = 1;
    nx = length(ctrlKite.sys.StateName);
    for ix = 1:nx
        subplot(nx, 1, ix);  hold on, grid on
        title(ctrlKite.sys.StateName{ix}); ylabel(ctrlKite.sys.StateUnit{ix});
        % Multiple shooting
        plot(single{1}.X{jRun}.Time, single{1}.X{jRun}.(ix), 'o-', 'MarkerSize', 6, 'LineWidth', 1, 'DisplayName', 'mshoot'),
        % Multiple shooting with cosine grid
%         plot(single{2}.Xcos{jRun}.Time, single{2}.Xcos{jRun}.(ix), 'o-', 'MarkerSize', 6, 'LineWidth', 1, 'DisplayName', 'ms\_cos')
        % Collocation
        % %         plot(single{3}.X{jRun}.Time, single{3}.X{jRun}.(ix), 'LineWidth', 2, 'DisplayName', 'colloc')
        % %         plot(single{3}.Xcol{jRun}.Time, single{3}.Xcol{jRun}.(ix), '^k', 'MarkerSize', 6, 'LineWidth', 1, 'DisplayName', 'c\_points');
        
        if exist('single_cpp')
            % Multiple shooting
            plot(single_cpp.X.Time, single_cpp.X.(ix), 'o-', 'MarkerSize', 6, 'LineWidth', 1, 'DisplayName', 'cpp mshoot'),
        end
    end
    subplot(nx, 1, 1); legend('-DynamicLegend')
    clear runToPlot ix nx jRun nRuns
    return
end

%% Simulate closed-loop system
% Choose controller for simulation
% % controller = controllers{1}; % Multiple shooting
% % controller = controllers{2}; % Multiple shooting with cosine grid
controller = controllers{3}; % Collocation

%clear controllers
show_open_loop = false;     % Compute only one trajectory

sim.dt = 0.005; % 200 Hz
sim.tf = 3;

% Construct simulator
simulator = Simulator(sim.dt, @simKite.dynamics);

% Create time vector and trajectories
sim.T      = 0:sim.dt:sim.tf;
sim.X      = zeros(simKite.nx, length(sim.T));
sim.X(:,1) = simKite.defaultState;
sim.U      = zeros(simKite.nu, length(sim.T)-1);
sim.Uid    = zeros(size(sim.T));

t_ctrl_next = 0;
iCtrl = 0;
complog = DurationLogger;

for iT = 1:length(sim.T)-1
    
    % Check if next control should be computed (at controller rate or after
    % computation delay of last computation). If chosen to show only one
    % open-loop trajectory, only compute once.
    if  sim.T(iT) >= t_ctrl_next && ~(show_open_loop && iT > 1)
        
        % Solve OCP and log computation time
        t0 = tic;
        controller.set_x0( sim.X(1:ctrlKite.nx, max(iT-1, 1)) );    % Most recent state before start of computation (Previous or first in simulation)
        controller = controller.solve();
        iCtrl = iCtrl + 1;                                          % Count OCP solution
        comp_time = toc(t0);
        complog = complog.log(comp_time);
        
        % Determine at which simulation time to apply the control
        t_apply = t_ctrl_next + ctrlset.comp_delay;                    % Apply time of this solution
        iT_apply = find(sim.T >= t_apply, 1);                       % Index of apply time in simulation time grid
        
        % Apply control to simulation
        if ~isempty(iT_apply)                                       % If this is empty, trajectory would be applied after simulation end
            % Option a) Hold first value of trajectory
            % %             u = controller.get_u0();
            % %             u = max( min(u, vec(simKite.phyUBU)), vec(simKite.phyLBU));       % Clip to physical inputs
            % %             sim.U(:, iT_apply:end) = repmat(u, 1, size(sim.U, 2)-iT_apply+1); % Insert first controller action until the end
            % %             sim.Uid(:, iT_apply:end) = iCtrl;                                 % Record where in grid which trajectory has been applied
            
            % Option b) Follow control trajectory during computation of next one
            Ui = controller.get_control_trajectory(sim.dt);                              % Get control trajectory resampled on simulation rate
            Ui = max( min(Ui, vec(simKite.phyUBU)), vec(simKite.phyLBU));                % Clip to physical inputs
            remaining_length = min( size(sim.U, 2)-iT_apply, size(Ui, 2) );              % Determine length to fill in simulation grid (end of trajectory or simulation)
            sim.U(:, iT_apply:iT_apply+remaining_length-1 ) = Ui(:, 1:remaining_length); % Insert full trajectory in remaining simulation control
            sim.Uid(:, iT_apply:iT_apply+remaining_length-1 ) = iCtrl;                   % Record where in grid which trajectory has been applied
            
            % Determine start time of next computation (scheduled controller dt or computation delay)
            t_ctrl_next = t_ctrl_next + max(ctrlset.dt, ctrlset.comp_delay);
        end
    end
    
    % Simulate system
    sim.X(:,iT+1) = simulator.simulate_step( sim.X(:,iT), sim.U(:,iT) );
    
    % Add noise
    % %     sim.X(:,iT+1) = sim.X(:,iT+1) + [0.01; deg2rad(0.2); deg2rad(0.5); deg2rad(0); 0; 0] .* randn(simKite.nx, 1); % Add noise
    
end, clear iT iCtrl remaining_length x u Ui t_apply iT_apply t_ctrl_next t0 comp_time

complog.printStatistics();

%% Plot trajectory
plot_hdls = simulator.plot_control_state_trajectory(sim.T, sim.U, sim.X, simKite.sys);

% Plot point at begin of each new trajectory
newUid = sim.Uid;
newUid( diff([0 sim.Uid]) == 0 ) = 0;
plot(plot_hdls.ax(1), sim.T(newUid ~= 0), rad2deg(sim.U(1, newUid ~= 0)), '.k', 'MarkerSize', 8); clear newUid

if ctrlobj.track_angle
    line(plot_hdls.ax(6), xlim(plot_hdls.ax(6)), rad2deg([ref.angle ref.angle]), 'LineWidth', 1, 'Color','black','LineStyle','--')
end
if ctrlobj.track_posDown
    line(plot_hdls.ax(8), xlim(plot_hdls.ax(8)), [ref.h ref.h], 'LineWidth', 1, 'LineWidth', 1, 'Color','black','LineStyle','--')
end
