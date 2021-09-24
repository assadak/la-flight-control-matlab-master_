classdef LonKiteMsOcp
    
    properties (Access = public)
        tf          % Final time of trajectory
        Ts          % Sample time
        N           % Number of points
        
        use_cos_grid % (Logical) Use cosine time grid
        T           % Time trajectory (times of grid points)
        
        opti        % Optimization object
        sol         % Solution object
    end
    
    properties (Access = public)
        control_objectives % Logical array encoding control objectives
        
        % Decision variables
        X_sym            % State trajectory in matrix shape (nx, N+1)
        U_sym            % Control trajectory in matrix shape (nu, N)
        
        x0_sym           % Initial state (nx, 1)
        
        slacks
        
        sym_p            % Symbolic parameter struct (setup in constructor)
    end
    
    methods (Access = public)
        function obj = LonKiteMsOcp(tf, N, lonKite, control_objectives, use_cos)
            
            % Default arguments
            if nargin < 5
                control_objectives = struct(...
                    'track_angle',      1, ...
                    'track_airspeed',   0, ...
                    'minimize_control', 1, ...
                    'track_posDown',    0);
            end
            if nargin < 4
                use_cos = false;
            end
            obj.control_objectives = control_objectives;
            
            % Save trajectory parameters
            obj.tf = tf;
            obj.N = N;
            obj.Ts = tf/N;
            
            fprintf(['Multiple shooting - N: ' num2str(N) ', Ts: ' num2str(obj.Ts) ', use cos: ' num2str(use_cos) '\n'])
            
            obj.use_cos_grid = use_cos;
            if use_cos
                obj.T = obj.tf * obj.get_cos_grid(obj.N+1);
            else
                obj.T = linspace(0, obj.tf, N+1);
            end
            %T = obj.T
            
            % Create casadi OptiStack
            import casadi.*
            opti = casadi.Opti(); % Optimization problem
            
            n = length(lonKite.sys.StateName);  % Number of states
            m = length(lonKite.sys.InputName);  % Number of inputs
            
            % ---- decision variables ---------
            X  = opti.variable(n, N+1); % State trajectory variables
            U  = opti.variable(m, N);   % Control trajectory variables
            x0 = opti.parameter(n, 1);  % Initial state
            
            obj.sym_p.config.angle_ref    = opti.parameter(1, 1);  % Reference angle
            obj.sym_p.config.Va_ref       = opti.parameter(1, 1);  % Reference airspeed
            obj.sym_p.config.h_ref        = opti.parameter(1, 1);  % Reference height
            
            obj.sym_p.tuning.mayer_multiplier = opti.parameter(1, 1);  % Multiplier for non_control stage cost to obtain mayer cost
            
            obj.sym_p.tuning.W_angle_err  = opti.parameter(1, 1);  % Weight on angle error
            obj.sym_p.tuning.W_Va_err     = opti.parameter(1, 1);  % Weight on airspeed error
            obj.sym_p.tuning.W_h_err      = opti.parameter(1, 1);  % Weight on height error
            
            obj.sym_p.tuning.R_diag       = opti.parameter(m, 1);  % Diagonal of control weight matrix
            
            
            % ---- objective ---------
            hs = diff(obj.T);  % Step durations
            cost = 0;
            for k = 1:N
                % Stage cost
                cost = cost + hs(k)/mean(hs) * (...
                    obj.get_non_control_cost( X(:,k) ) + ...           % Standard formulation
                    ... k/N * obj.get_non_control_cost( X(:,k) ) + ... % Later deviation from objective is more expensive: Act earlier!
                    obj.get_control_cost( U(:,k) ) );
            end
            
            % Terminal cost (Meyer term)
            cost = cost + obj.sym_p.tuning.mayer_multiplier * obj.get_non_control_cost( X(:,N+1)  );
            
            opti.minimize( cost );
            
            % ---- dynamics --------
            %             f_discrete = @(x_,u_,p_) RK2(x_, u_, p_, Ts, @lonKite.dynamics);
            f_discrete = @(x_,u_,p_,h_) RK4(x_, u_, p_, h_, @lonKite.dynamics);
            for k = 1:N
                opti.subject_to( X(:,k+1) == f_discrete(X(:,k), U(:,k), 0, hs(k)) );
            end
            
            % ---- boundary conditions --------
            opti.subject_to( X(:,1) == x0 );   % use initial position
            
            % ---- Setup solver NLP    ------
            opts = struct('ipopt', struct('print_level', 0, 'tol', 1e-3), 'print_time', false);
            opti.solver('ipopt', opts);
            
            
            % Save variables to object
            obj.opti = opti;
            obj.X_sym = X;
            obj.U_sym = U;
            obj.x0_sym = x0;
        end
        
        % Setters for config parameters
        function set_angle_ref(obj, angle_ref)
            obj.opti.set_value( obj.sym_p.config.angle_ref, angle_ref);
        end
        function set_Va_ref(obj, Va_ref)
            obj.opti.set_value( obj.sym_p.config.Va_ref, Va_ref);
        end
        function set_h_ref(obj, h_ref)
            obj.opti.set_value( obj.sym_p.config.h_ref, h_ref);
        end
        
        % Setters for tuning parameters
        function set_mayer_multiplier(obj, mayer_multiplier)
            obj.opti.set_value( obj.sym_p.tuning.mayer_multiplier, mayer_multiplier);
        end
        
        function set_W_Va_err(obj, W_Va_err)
            obj.opti.set_value( obj.sym_p.tuning.W_Va_err, W_Va_err);
        end
        function set_W_angle_err(obj, W_angle_err)
            obj.opti.set_value( obj.sym_p.tuning.W_angle_err, W_angle_err);
        end
        function set_W_h_err(obj, W_h_err)
            obj.opti.set_value( obj.sym_p.tuning.W_h_err, W_h_err);
        end
        
        function set_R_diag(obj, R_diag)
            obj.opti.set_value( obj.sym_p.tuning.R_diag, R_diag);
        end
        
        % Setters for state and control bounds
        function obj = set_ubx(obj, ubx, slack_bounds, slack_weights)
            
            ubx = [vec(ubx); inf( size(obj.X_sym,1)-length(ubx), 1 )];
            
            if nargin < 4
                if nargin < 3
                    slack_bounds = zeros(size(ubx));
                end
                slack_weights = 10000 * ones(size(slack_bounds));
            end
            slack_bounds = [vec(slack_bounds); zeros( size(obj.X_sym,1)-length(slack_bounds), 1 )];
            slack_weights = [vec(slack_weights); zeros( size(obj.X_sym,1)-length(slack_weights), 1 )];
            
            for ix = 1:length(ubx)
                
                if (ubx(ix) < inf && slack_bounds(ix) > 0)
                    % Create slack vector variable
                    slack = obj.opti.variable(1, 1);
                    
                    % Set constraint
                    obj.opti.subject_to( obj.X_sym(ix,:) <= ubx(ix) + slack );
                    obj.opti.subject_to( 0 <= slack <= slack_bounds(ix) );
                    
                    % Add penalty to cost function
                    cost = obj.opti.f();
                    slack_penalty = slack_weights(ix) * 0.5*(slack^2 + slack);
                    obj.opti.minimize( cost + slack_penalty );
                    
                    obj.slacks = [obj.slacks; slack];
                    disp(['Set soft upper constraint for state no ' num2str(ix) ': ' num2str(ubx(ix)) ' (+ ' num2str(slack_bounds(ix)) ')']);
                else
                    % Hard constraints
                    disp(['Set hard upper constraint for state no ' num2str(ix) ': ' num2str(ubx(ix))]);
                    obj.opti.subject_to( obj.X_sym(ix,:) <= ubx(ix) );
                end
            end
        end
        function obj = set_lbx(obj, lbx, slack_bounds, slack_weights)
            
            lbx = [vec(lbx); -inf( size(obj.X_sym,1)-length(lbx), 1 )];
            
            if nargin < 4
                if nargin < 3
                    slack_bounds = zeros(size(lbx));
                end
                slack_weights = 10000 * ones(size(slack_bounds));
            end
            slack_bounds = [vec(slack_bounds); zeros( size(obj.X_sym,1)-length(slack_bounds), 1 )];
            slack_weights = [vec(slack_weights); zeros( size(obj.X_sym,1)-length(slack_weights), 1 )];
            
            for ix = 1:length(lbx)
                
                if (lbx(ix) < inf && slack_bounds(ix) > 0)
                    % Create slack vector variable
                    slack = obj.opti.variable(1, 1);
                    
                    % Set constraint
                    obj.opti.subject_to( lbx(ix) - slack <= obj.X_sym(ix,:) );
                    obj.opti.subject_to( 0 <= slack <= slack_bounds(ix) );
                    
                    % Add penalty to cost function
                    cost = obj.opti.f();
                    slack_penalty = slack_weights(ix) * 0.5*(slack^2 + slack);
                    obj.opti.minimize( cost + slack_penalty );
                    
                    obj.slacks = [obj.slacks; slack];
                    disp(['Set soft lower constraint for state no ' num2str(ix) ': ' num2str(lbx(ix)) ' (- ' num2str(slack_bounds(ix)) ')']);
                else
                    % Hard constraints
                    disp(['Set hard lower constraint for state no ' num2str(ix) ': ' num2str(lbx(ix))]);
                    obj.opti.subject_to( lbx(ix) <= obj.X_sym(ix,:) );
                end
            end
        end
        function set_ubu(obj, ubu)
            obj.opti.subject_to( obj.U_sym(:) <= repmat( vec(ubu), obj.N, 1) );
        end
        function set_lbu(obj, lbu)
            obj.opti.subject_to( repmat( vec(lbu), obj.N, 1) <= obj.U_sym(:) );
        end
        
        % Setters for initial guesses
        function set_X_guess(obj, X_guess)
            obj.opti.set_initial(obj.X_sym, X_guess);
        end
        function set_U_guess(obj, U_guess)
            obj.opti.set_initial(obj.U_sym, U_guess);
        end
        
        % Setter for initial state
        function set_x0(obj, x0)
            obj.opti.set_value( obj.x0_sym, x0 );
        end
        
        % Call for OCP solve
        function [obj] = solve(obj)
            
            if nargin < 2
                warm_start = true;
            end
            
            obj.sol = obj.opti.solve();
            assert(obj.sol.stats.success == 1, 'Error computing optimal input');
            
            if warm_start
                % Use the current solution to speed up the next optimization
                obj.opti.set_initial( obj.sol.value_variables() );
                obj.opti.set_initial( obj.opti.lam_g, obj.sol.value(obj.opti.lam_g) );
            end
        end
        
        % Getters for solution parts
        function u = get_u0(obj)
            u = obj.sol.value( obj.U_sym(:,1) );
        end
        function [U, T, Uopt, Topt] = get_control_trajectory(obj, dt)
            
            Topt = obj.T(1:end-1);
            Uopt = obj.sol.value( obj.U_sym );
            [U, T] = obj.wrap_retime(Uopt, Topt, dt);
        end
        function [X, T, Xopt, Topt] = get_state_trajectory(obj, dt)

            Topt = obj.T;
            Xopt = obj.sol.value( obj.X_sym );
            [X, T] = obj.wrap_retime(Xopt, Topt, dt);
        end
    end
    
    methods (Access = private)
        % Helpers for non-/control stage costs
        function non_control_cost = get_non_control_cost(obj, x)
            
            non_control_cost = 0;
            
            if obj.control_objectives.track_angle
                % Objective: Track angle
                angle_err = obj.sym_p.config.angle_ref - x(4);
                non_control_cost = non_control_cost + obj.sym_p.tuning.W_angle_err * angle_err^2;
            end
            
            if obj.control_objectives.track_airspeed
                % Objective: Track airspeed
                Va_err = obj.sym_p.config.Va_ref - x(1);
                non_control_cost = non_control_cost + obj.sym_p.tuning.W_Va_err * Va_err^2;
            end
            
            if obj.control_objectives.track_posDown
                height_err = obj.sym_p.config.h_ref -(-x(6));
                non_control_cost = non_control_cost + obj.sym_p.tuning.W_h_err * height_err^2;
            end
        end
        
        function control_cost = get_control_cost(obj, u)
            
            control_cost = 0;
            
            if obj.control_objectives.minimize_control
                % Objective: Minimize control effort
                control_cost = control_cost + u' * diag(obj.sym_p.tuning.R_diag) * u;
            end
        end
    end
    
    methods (Static)
        function Tcos = get_cos_grid(N)
            Tcos = (1 - cos( pi/2 * linspace(0, 1, N)));
        end
        
        function [X, T] = wrap_retime(X, Tx, dt)
            Xtt = timetable(seconds(Tx)', X');                           
            Xrt = retime(Xtt, 'regular', 'previous', 'TimeStep', seconds(dt));
            X = Xrt.Variables';
            T = seconds(Xrt.Time)';
        end
    end
end

