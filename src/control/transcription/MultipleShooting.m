classdef MultipleShooting < TranscriptionBase
    %MULTIPLESHOOTING Summary of this class goes here
    %   Detailed explanation goes here
    
    methods (Access = public)
        % Override abstract base class methods
        function [ocp, obj] = generate_decision_vars(obj)
            
            ocp = obj.ocp;
            nx = ocp.model.nx;  % Number of states
            nu = ocp.model.nu;  % Number of inputs
            
            N = obj.opts.N;
            ocp.N = N;
            ocp.Nseg = N;
            
            % Generate trajectory state time grid
            if obj.opts.use_cos_grid
                ocp.Tx = ocp.tf * obj.get_cos_grid(ocp.N+1);
            else
                ocp.Tx = linspace(0, ocp.tf, N+1);
            end
            ocp.Tu = ocp.Tx(1:end-1);

            % Create symbolic decision variables
            ocp.X_sym = ocp.opti.variable(nx, N+1); % State trajectory variables
            ocp.U_sym = ocp.opti.variable(nu, N); % Control trajectory variables
            ocp.x0_sym = ocp.opti.parameter(nx, 1); % Initial state
        end
        
        function setup_cost(obj)
            
            N = obj.ocp.N;
            X = obj.ocp.X_sym;
            U = obj.ocp.U_sym;
            
            hs = diff(obj.ocp.Tx);  % Step durations

            cost = 0;
            for i = 0:N-1
                cost = cost + hs(i+1)/mean(hs) * obj.ocp.stage_cost_impl( X(:,i+1), U(:,i+1) );
            end
            
            cost = cost + obj.ocp.terminal_cost_impl( X(:,N+1) );
            
            obj.ocp.opti.minimize( cost );
        end
        
        function setup_dyn_constraints(obj)
            
            ocp = obj.ocp;
            model = ocp.model;
            f_discrete = @(x_,u_,p_,h_) RK4(x_, u_, p_, h_, @model.dynamics);
            N = ocp.N;
            X = ocp.X_sym;
            U = ocp.U_sym;
            hs = diff(ocp.Tx);
            
            for i = 0:N-1
                obj.ocp.opti.subject_to( X(:,(i+1)+1) == f_discrete(X(:,i+1), U(:,i+1), 0, hs(i+1)) );
            end
            
        end
        
        % Output functions
        % Use base class implementations / do not override
    end
    
    methods (Static, Access = private)
        function Tcos = get_cos_grid(N)
            Tcos = (1 - cos( pi/2 * linspace(0, 1, N)));
        end
    end
end
