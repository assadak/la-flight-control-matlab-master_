classdef LonKiteColOcp < Approximation
    
    properties (Access = public)
        N           % Total number of points (Nseg * Dpoly + 1)
        T           % Times at collocation points
        N_nodes_per_segment
        
        control_objectives % Logical array encoding control objectives
        Nu                 % Control trajectory points
        approx_control     % Logical: Polynomial approximation of continuous control
                
        opti        % Optimization object
        sol         % Solution object
    end
    
    properties (Access = protected)
        % Decision variables
        X_sym            % State trajectory in matrix shape (nx, N+1)
        U_sym            % Control trajectory in matrix shape (nu, N)
        
        x0_sym           % Initial state (nx, 1)
        
        slacks
        
        sym_p            % Symbolic parameter struct (setup in constructor)
    end
    
    methods (Access = public)
        function [obj, N] = LonKiteColOcp(tf, Nseg, D, lonKite, control_objectives, approx_control)

            if D < 1
                fprintf('Polynomial degree must be > 0. Continuing with D = 1.\n')
                D = 1;
            end
            obj = obj@Approximation(tf, Nseg, D);
            d_seg = obj.d_seg;
            Tau = obj.Tau;
            
            fprintf(['Collocation       - H: ' num2str(tf) ', Nseg: ' num2str(Nseg) ', Dpoly: ' num2str(D) '\n'])
            
            % Default arguments
            if nargin < 6
                approx_control = false;
            end
            obj.approx_control = approx_control;
 
            if nargin < 5
                control_objectives = struct(...
                    'track_angle',      1, ...
                    'track_airspeed',   0, ...
                    'minimize_control', 1, ...
                    'track_posDown',    0);
            end
            obj.control_objectives = control_objectives;
            
            % Determine number of state trajectory points
            N_nodes_per_segment = D+1; % 0, ..., D
            N = Nseg * N_nodes_per_segment + 1;
            %(0:N-1)
            
            % Generate trajectory state time grid
            T = 0:N-1; % For debugging
            for iSeg = 0:Nseg-1
                % Segment indices (0-based)
                i0 = iSeg * N_nodes_per_segment;
                iD = (iSeg+1) * N_nodes_per_segment - 1;
                T((i0:iD)+1) = iSeg * d_seg + d_seg * Tau;
            end
            T(end) = tf;
            
            % Determine number of control trajectory points
            if approx_control
                Nu = N;
            else
                Nu = Nseg;
            end
            
            % Save to object
            obj.N_nodes_per_segment = N_nodes_per_segment;
            obj.N = N;
            obj.Nu = Nu;
            obj.T = T;
            
            % ---- Formulate OCP ------------------------------------------
            % Create casadi OptiStack
            import casadi.*
            opti = casadi.Opti(); % Optimization problem
            obj.opti = opti;
            
            nx = length(lonKite.sys.StateName);  % Number of states
            nu = length(lonKite.sys.InputName);  % Number of inputs
            
            % ---- decision variables -------------------------------------
            X = opti.variable(nx, N); % State trajectory variables
            U = opti.variable(nu, Nu); % Control trajectory variables
            x0 = opti.parameter(nx, 1); % Initial state
            obj.X_sym = X;
            obj.U_sym = U;
            obj.x0_sym = x0;
            
            obj.sym_p.config.angle_ref    = opti.parameter(1, 1);  % Reference angle
            obj.sym_p.config.Va_ref       = opti.parameter(1, 1);  % Reference airspeed
            obj.sym_p.config.h_ref        = opti.parameter(1, 1);  % Reference height
            
            obj.sym_p.tuning.mayer_multiplier = opti.parameter(1, 1);  % Multiplier for non_control stage cost to obtain mayer cost
            
            obj.sym_p.tuning.W_angle_err  = opti.parameter(1, 1);  % Weight on angle error
            obj.sym_p.tuning.W_Va_err     = opti.parameter(1, 1);  % Weight on airspeed error
            obj.sym_p.tuning.W_h_err      = opti.parameter(1, 1);  % Weight on height error
            
            obj.sym_p.tuning.R_diag       = opti.parameter(nu, 1);  % Diagonal of control weight matrix
            
            % ---- objective ----------------------------------------------
            % >>>>>> For comparison only!
% %             disp('Sumarize stage costs at collocation points') % !For comparison only!
% %             cost = 0;
% %             for i = 0:N-2
% %                 cost = cost + obj.get_non_control_cost( X(:,i+1) ) + obj.get_control_cost( U(:,i+1) );
% %             end
            % <<<<<< For comparison only!
            
            disp('Add integral cost')
            cost = obj.get_integral_cost(X, U);
            
            % Terminal cost (Meyer term)
            disp('Add terminal cost')
            cost = cost + obj.sym_p.tuning.mayer_multiplier * obj.get_non_control_cost( X(:,(N-1)+1)  );
            
            opti.minimize( cost );
            
            % ---- dynamics -----------------------------------------------
            fprintf('\nSet dynamics and multi-segment constraints:\n')
            
            disp('seg 0, point 0 (traj point 0): x = x0')
            opti.subject_to( X(:,1) == x0 );  % X(:,0) = x0

            for iSeg = 0:Nseg-1
                
                Xseg = obj.get_segment(X, iSeg);
                if obj.approx_control
                    Useg = obj.get_segment(U, iSeg); % Get collocation values of this segment
                else
                    Useg = repmat( U(:,iSeg+1), 1, size(Xseg,2) ); % Duplicate single segment value for same interfacing
                end
                
                % Enforce dynamics at each collocation point
                for jD = 1:D
                    i = iSeg * (N_nodes_per_segment) + jD; disp(['seg ' num2str(iSeg) ', point ' num2str(jD) ' (traj point ' num2str(i) '): dx = f(x, u)'])
                    opti.subject_to( obj.dapprox(Xseg, Tau(jD+1)) == lonKite.dynamics( Xseg(:,jD+1), Useg(:,jD+1), 0) );
                end
                
                if iSeg < Nseg-1
                    % Glue segments together
                    fprintf(['--- seg border ' num2str(iSeg) '|' num2str(iSeg+1) ':\n']);
                    
                    Xseg_next = obj.get_segment(X, iSeg+1);
                    opti.subject_to( Xseg_next(:,0+1) == obj.approx(Xseg, 1) );     fprintf(['       x_' num2str(iSeg+1) ',0 = approx(Xseg_' num2str(iSeg) ', 1)\n']);
                    
                    if approx_control
                        Useg_next = obj.get_segment(U, iSeg+1);
                        opti.subject_to( Useg_next(:,0+1) == obj.approx(Useg, 1) ); fprintf(['       u_' num2str(iSeg+1) ',0 = approx(Useg_' num2str(iSeg) ', 1)\n']);
                        opti.subject_to( obj.dapprox(Useg_next, 0) == obj.dapprox(Useg, 1) ); fprintf(['       approx(Useg_' num2str(iSeg+1) ', 0) = approx(Useg_' num2str(iSeg) ', 1)\n']);
                    end
                    fprintf(['---------------------------------------------\n'])
                else
                    fprintf(['--- seg border ' num2str(iSeg) ':\n']);
                    
                    opti.subject_to( X(:,end) == obj.approx(Xseg, 1) );             fprintf(['       x_' num2str(iSeg+1) ',0 = approx(Xseg_' num2str(iSeg) ', 1)\n']);
                    
                    if approx_control
                        opti.subject_to( U(:,end) == obj.approx(Useg, 1) );         fprintf(['       u_' num2str(iSeg+1) ',0 = approx(Useg_' num2str(iSeg) ', 1)\n']);
%                         opti.subject_to( obj.dapprox(obj.get_segment(U, Nseg-1), 1) == 0 ); fprintf(['       approx(Useg_' num2str(iSeg) ', 1) = 0\n']);
                    end
                    fprintf(['---------------------------------------------\n'])
                end
            end
            
            % ---- Setup solver NLP    ------
            opts = struct('ipopt', struct('print_level', 0, 'tol', 1e-3), 'print_time', false);
            opti.solver('ipopt', opts);
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
            obj.opti.subject_to( obj.U_sym(:) <= repmat( vec(ubu), obj.Nu, 1) );
        end
        function set_lbu(obj, lbu)
            obj.opti.subject_to( repmat( vec(lbu), obj.Nu, 1) <= obj.U_sym(:) );
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
            
            Uopt = obj.sol.value( obj.U_sym );
            
            if obj.approx_control
                % Control is stored as collocation points -> resample
                Topt = obj.T;
                [U, T] = obj.resample_trajectory(Uopt, dt);
            else
                % Control is piecewise constant on segments
                Topt = obj.Tsegs;
                [U, T] = obj.wrap_retime(Uopt, Topt, dt);
            end
        end
        function [X, T, Xopt, Topt] = get_state_trajectory(obj, dt)
            
            Xopt = obj.sol.value( obj.X_sym );
            Topt = obj.T;
            
            % Resample whole trajectory
            [X, T] = obj.resample_trajectory(Xopt, dt);
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
        
        % Compute integral cost
        function cost = get_integral_cost(obj, X, U)
            
            Nseg = obj.Nseg;
            D = obj.D;
            Tau = obj.Tau;
            
            colCostVectSum = zeros(D, 1);
            
            for iSeg = 0:Nseg-1
                
                Xseg = obj.get_segment(X, iSeg);
                if obj.approx_control
                    Useg = obj.get_segment(U, iSeg); % Get collocation values of this segment
                else
                    Useg = repmat( U(:,iSeg+1), 1, size(Xseg,2) ); % Duplicate single segment value for same interfacing
                end
                
                % Build vector where rows contain cost at jD = 1 ... D collocation points
                colCostVect = casadi.MX(D, 1);
                for jD = 1:D
                    colCostVect(jD) = obj.get_non_control_cost( Xseg(:,jD+1) ) + obj.get_control_cost( Useg(:,jD+1) );
                end
                
                colCostVectSum = colCostVectSum + colCostVect;
            end
            
            L_1D_1 = zeros(1, D);
            for i = 1:D
                L_1D_1(i) = obj.L(i, 1);
            end
            
            dL_mat = zeros(D, D);
            for i = 1:D
                for j = 1:D
                    dL_mat(i, j) = obj.dL(j, Tau(i+1));
                end
            end
            
            LAMBDA_T = L_1D_1 * dL_mat^-1;
            
            cost = obj.d_seg * LAMBDA_T * colCostVectSum;
        end
        
        
        % Data handling
        function SEG = get_segment(obj, TRAJ, iSeg)
            % Assume trajectory grid is like illustrated below:
            % (see http://casadi.sourceforge.net/v1.9.0/users_guide/html/node7.html#SECTION00740000000000000000)
            % traj point    0  1  2 ...                                            nSeg*(D+1)+1
            % seg  point    0  1  ...  D                0  1  ...  D
            %                             0  1  ...  D                0  1  ...  D
            %               |             |             |             |             |
            % iSeg          0             1            ...                         Nseg-1
            %
            %             Return ith SEG |____________|
            
            seg0 = iSeg * obj.N_nodes_per_segment;
            segE = (iSeg+1) * obj.N_nodes_per_segment - 1;
            
            SEG = TRAJ(:, (seg0:segE)+1);
        end
        
        function [X, Tn] = resample_trajectory(obj, Xcoll, dt)
            
            Tn = 0:dt:obj.tf; % New, resampled time grid
            X = zeros(size(Xcoll, 1), length(Tn));
            t_seg = -2*obj.d_seg;
            iSeg = -1;
            
            for i = 1:length(Tn)
                if Tn(i) > (t_seg + obj.d_seg + 0.00001)
                    % Set next segment for resampling
                    iSeg = iSeg + 1;
                    t_seg = iSeg * obj.d_seg;
                    Xseg = obj.get_segment(Xcoll, iSeg);
                end
                % Resample from this segment
                X(:, i) = obj.approx(Xseg, (Tn(i)-t_seg)/obj.d_seg);
            end
        end
        
        function [X, T] = wrap_retime(obj, X, Tx, dt)
            Xtt = timetable(seconds(Tx)', X');
            Xrt = retime(Xtt, 'regular', 'previous', 'TimeStep', seconds(dt));
            X = Xrt.Variables';
            T = seconds(Xrt.Time)';
        end
    end
    
end

