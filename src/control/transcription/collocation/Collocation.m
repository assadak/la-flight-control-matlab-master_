classdef Collocation < TranscriptionBase
    %COLLOCATION Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access = public)
        approx
        N_nodes_per_segment
    end
    
    methods (Access = public)
        % Override abstract base class methods
        function [ocp, obj] = generate_decision_vars(obj)

            ocp = obj.ocp;
            opts = obj.opts;
            
            nx = ocp.model.nx;  % Number of states
            nu = ocp.model.nu;  % Number of inputs
            tf = ocp.tf;
            
            Nseg = opts.Nseg;
            D    = opts.D;
            approx_control = opts.approx_control;
            
            if D < 1
                fprintf('Polynomial degree must be > 0. Continuing with D = 1.\n')
                D = 1;
            end
            obj.approx = PolyApprox(ocp.tf, Nseg, D);
            
            fprintf(['Collocation       - H: ' num2str(tf) ', Nseg: ' num2str(Nseg) ', Dpoly: ' num2str(D) '\n'])
            
            % Determine number of state trajectory points
            obj.N_nodes_per_segment = D+1; % 0, ..., D
            N_nodes_per_segment =  obj.N_nodes_per_segment;
            N = Nseg * N_nodes_per_segment + 1;
            
            ocp.Nseg = Nseg;
            ocp.N = N;
            %(0:N-1)
            
            d_seg = obj.approx.d_seg;
            Tau = obj.approx.Tau;
            
            % Generate trajectory state time grid
            Tx = 0:N-1; % For debugging
            for iSeg = 0:Nseg-1
                % Segment indices (0-based)
                i0 = iSeg * N_nodes_per_segment;
                iD = (iSeg+1) * N_nodes_per_segment - 1;
                Tx((i0:iD)+1) = iSeg * d_seg + d_seg * Tau;
            end
            Tx(end) = tf;
            
            if approx_control
                Nu = N;
                Tu = Tx;
            else
                Nu = Nseg;
                Tu = (0:Nseg-1)*d_seg;
            end
            ocp.Tx = Tx;
            ocp.Tu = Tu;
            
            % Create symbolic decision variables
            ocp.X_sym = ocp.opti.variable(nx, N); % State trajectory variables
            ocp.U_sym = ocp.opti.variable(nu, Nu); % Control trajectory variables
            ocp.x0_sym = ocp.opti.parameter(nx, 1); % Initial state
        end
        
        function setup_cost(obj)
            
            N = obj.ocp.N;
            X = obj.ocp.X_sym;
            U = obj.ocp.U_sym;
            
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
            cost = cost + obj.ocp.terminal_cost_impl( X(:,(N-1)+1) );
            
            obj.ocp.opti.minimize( cost );
        end
        
        function setup_dyn_constraints(obj)
            
            ocp = obj.ocp;
            model = ocp.model;
            opts = obj.opts;
            
            X = ocp.X_sym;
            U = ocp.U_sym;
            
            apr = obj.approx;
            Nseg = apr.Nseg;
            D    = apr.D;
            Tau = apr.Tau;
            
            approx_control = opts.approx_control;
            
            fprintf('\nSet dynamics and multi-segment constraints:\n')
            
            %             disp('seg 0, point 0 (traj point 0): x = x0')
            %             opti.subject_to( X(:,1) == x0 );  % X(:,0) = x0
            
            for iSeg = 0:Nseg-1
                
                Xseg = obj.get_segment(X, iSeg);
                if approx_control
                    Useg = obj.get_segment(U, iSeg); % Get collocation values of this segment
                else
                    Useg = repmat( U(:,iSeg+1), 1, size(Xseg,2) ); % Duplicate single segment value for same interfacing
                end
                
                % Enforce dynamics at each collocation point
                for jD = 1:D
                    i = iSeg * obj.N_nodes_per_segment + jD; disp(['seg ' num2str(iSeg) ', point ' num2str(jD) ' (traj point ' num2str(i) '): dx = f(x, u)'])
                    ocp.opti.subject_to( apr.dapprox(Xseg, Tau(jD+1)) == model.dynamics( Xseg(:,jD+1), Useg(:,jD+1), 0) );
                end
                
                if iSeg < Nseg-1
                    % Glue segments together
                    fprintf(['--- seg border ' num2str(iSeg) '|' num2str(iSeg+1) ':\n']);
                    
                    Xseg_next = obj.get_segment(X, iSeg+1);
                    ocp.opti.subject_to( Xseg_next(:,0+1) == apr.approx(Xseg, 1) );     fprintf(['       x_' num2str(iSeg+1) ',0 = approx(Xseg_' num2str(iSeg) ', 1)\n']);
                    
                    if approx_control
                        Useg_next = obj.get_segment(U, iSeg+1);
                        ocp.opti.subject_to( Useg_next(:,0+1) == apr.approx(Useg, 1) ); fprintf(['       u_' num2str(iSeg+1) ',0 = approx(Useg_' num2str(iSeg) ', 1)\n']);
                        ocp.opti.subject_to( apr.dapprox(Useg_next, 0) == apr.dapprox(Useg, 1) ); fprintf(['       approx(Useg_' num2str(iSeg+1) ', 0) = approx(Useg_' num2str(iSeg) ', 1)\n']);
                    end
                    fprintf(['---------------------------------------------\n'])
                else
                    fprintf(['--- seg border ' num2str(iSeg) ':\n']);
                    
                    ocp.opti.subject_to( X(:,end) == apr.approx(Xseg, 1) );             fprintf(['       x_' num2str(iSeg+1) ',0 = approx(Xseg_' num2str(iSeg) ', 1)\n']);
                    
                    if approx_control
                        ocp.opti.subject_to( U(:,end) == apr.approx(Useg, 1) );         fprintf(['       u_' num2str(iSeg+1) ',0 = approx(Useg_' num2str(iSeg) ', 1)\n']);
                        %                         opti.subject_to( obj.dapprox(obj.get_segment(U, Nseg-1), 1) == 0 ); fprintf(['       approx(Useg_' num2str(iSeg) ', 1) = 0\n']);
                    end
                    fprintf(['---------------------------------------------\n'])
                end
            end
            
        end
        
        % Override output functions with collocation-specific
        % implementations
        function [U, T] = get_control_trajectory_impl(obj, Uopt, Topt, dt)
            
            if obj.opts.approx_control
                % Control is stored on uneven time grid -> resample
                [U, T] = obj.resample_trajectory(Uopt, dt);
            else
                % Control is piecewise constant on segments
                [U, T] = obj.wrap_retime(Uopt, Topt, dt);
            end
        end
        function [X, T] = get_state_trajectory_impl(obj, Xopt, Z, dt)
            
            % Control is piecewise constant on segments
            [X, T] = obj.resample_trajectory(Xopt, dt);
        end
    end
    
    methods (Access = private)
        % Compute integral cost
        function cost = get_integral_cost(obj, X, U)
            
            opts = obj.opts;
            
            Nseg = opts.Nseg;
            D    = opts.D;
            Tau  = obj.approx.Tau;
            
            colCostVectSum = zeros(D, 1);
            
            for iSeg = 0:Nseg-1
                
                Xseg = obj.get_segment(X, iSeg);
                if obj.opts.approx_control
                    Useg = obj.get_segment(U, iSeg); % Get collocation values of this segment
                else
                    Useg = repmat( U(:,iSeg+1), 1, size(Xseg,2) ); % Duplicate single segment value for same interfacing
                end
                
                % Build vector where rows contain cost at jD = 1 ... D collocation points
                colCostVect = casadi.MX(D, 1);
                for jD = 1:D
                    colCostVect(jD) = obj.ocp.stage_cost_impl( Xseg(:,jD+1), Useg(:,jD+1) );
                end
                
                colCostVectSum = colCostVectSum + colCostVect;
            end
            
            L_1D_1 = zeros(1, D);
            for i = 1:D
                L_1D_1(i) = obj.approx.L(i, 1);
            end
            
            dL_mat = zeros(D, D);
            for i = 1:D
                for j = 1:D
                    dL_mat(i, j) = obj.approx.dL(j, Tau(i+1));
                end
            end
            
            LAMBDA_T = L_1D_1 * dL_mat^-1;
            
            cost = obj.approx.d_seg * LAMBDA_T * colCostVectSum;
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
            
            apr = obj.approx;
            d_seg = apr.d_seg;
            
            Tn = 0:dt:apr.tf; % New, resampled time grid
            X = zeros(size(Xcoll, 1), length(Tn));
            t_seg = -2 * d_seg;
            iSeg = -1;
            
            for i = 1:length(Tn)
                if Tn(i) > (t_seg + d_seg + 0.00001)
                    % Set next segment for resampling
                    iSeg = iSeg + 1;
                    t_seg = iSeg * d_seg;
                    Xseg = obj.get_segment(Xcoll, iSeg);
                end
                % Resample from this segment
                X(:, i) = apr.approx(Xseg, (Tn(i)-t_seg)/d_seg);
            end
        end
    end
    
end

