classdef OcpBase
    
    properties (Access = public)
        opts
        
        model
        nx
        nu
        
        tf          % End time
        N           % Total number of state trajectory points
        Tx          % Times at state trajectory points
        Tu          % Times at input trajectory points
        Nseg
        
        transcription
        opti        % Optimization object
        sol         % Solution object
    end
    properties (Access = public)
        % Decision variables
        X_sym            % State trajectory in matrix shape
        U_sym            % Control trajectory in matrix shape
        x0_sym           % Initial state
        p_sym            % Symbolic parameter struct
        
        slacks_sym
    end
    
    methods (Access = public)
        function obj = OcpBase(model, tf)
            obj.model = model;
            obj.nx = model.nx;
            obj.nu = model.nu;
            obj.tf = tf;
            
            import casadi.*
            opti = casadi.Opti(); % Optimization problem
            obj.opti = opti;
        end
        
        function obj = setup(obj, transcript_method)
            
            % Setup transcription
            if strcmp(transcript_method, 'collocation')
                transcr = Collocation();
                transcription_str = 'Collocation';
            else
                % Multiple shooting
                transcr = MultipleShooting();
                transcription_str = 'Multiple Shooting';
            end
            
            disp(['Transcribe OCP to NLP using ' transcription_str ' ...'])
            [obj, transcr] = transcr.transcribe(obj, obj.opts);
            obj.transcription = transcr;
            
            % Setup solver
            solver_opts = struct('ipopt', struct('print_level', 0, 'tol', 1e-3), 'print_time', false);
            obj.opti.solver('ipopt', solver_opts);
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
            obj.opti.subject_to( obj.U_sym(:) <= repmat(vec(ubu), size(obj.U_sym, 2), 1) );
        end
        function set_lbu(obj, lbu)
            obj.opti.subject_to( repmat(vec(lbu), size(obj.U_sym, 2), 1) <= obj.U_sym(:) );
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
            
            Topt = obj.Tu;
            Uopt = obj.sol.value( obj.U_sym );
            
            [U, T] = obj.transcription.get_control_trajectory_impl(Uopt, Topt, dt);
        end
        function [X, T, Xopt, Topt] = get_state_trajectory(obj, dt)
            
            Xopt = obj.sol.value( obj.X_sym );
            Topt = obj.Tx;
            
            [X, T] = obj.transcription.get_state_trajectory_impl(Xopt, Topt, dt);
        end
    end
    
    methods (Access = private)
        function cost = stage_cost_impl(obj, x, u)
            disp('Base class/cost_impl')
            cost = 0;
        end
        function cost = terminal_cost_impl(obj, x, u)
            disp('Base class/terminal_cost_impl')
            cost = 0;
        end
        function x_dot = cont_dynamics_impl(obj, x, u)
            disp('Base class/cont_dynamics_impl')
            x_dot = x;
        end
    end
    
end

