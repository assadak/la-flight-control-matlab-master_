classdef LonKiteDynamics < LonKiteParams
    
    properties (Access = public)
        
        nx = 4;         % Number of states
        nu = 2;         % Number of inputs
        stateRep = 'Pitch';
        
        defaultState    % State for which the dynamics is well defined
        defaultControl  % Control for which the dynamics is well defined
        phyUBU          % Physical upper control bound
        phyLBU          % Physical Lower control bound
        
        sys             % Struct containing system properties
     
        
    end
    
    % Public methods
    methods (Access = public)
        
        function obj = LonKiteDynamics(state_rep, model_params_filepath)
            
            obj.stateRep = state_rep;
            
            switch (obj.stateRep) %% chooses the state representation that is specified out of the 3 given
                
                case 'Flightpath'
                    sys.StateName = {'Va' 'alpha' 'q' 'gamma'};
                    sys.StateUnit = {'m/s' 'rad' 'rad/s' 'rad'};
                    sys.InputName = {'dE' 'dF'};
                    sys.InputUnit = {'rad' '-'};
                    
                case 'Pitch'
                    sys.StateName = {'Va' 'alpha' 'q' 'theta'};
                    sys.StateUnit = {'m/s' 'rad' 'rad/s' 'rad'};
                    sys.InputName = {'dE' 'dF'};
                    sys.InputUnit = {'rad' '-'};
                    
                case 'UW'
                    sys.StateName = {'vx' 'vz' 'q' 'theta'};
                    sys.StateUnit = {'m/s' 'm/s' 'rad/s' 'rad'};
                    sys.InputName = {'dE' 'dF'};
                    sys.InputUnit = {'rad' '-'};
                    
                otherwise
                    error(['State representation ''' obj.stateRep ''' not found.']);
            end
            obj.stateRep = state_rep;
            
            sys.OutputName = sys.StateName;
            sys.OutputUnit = sys.StateUnit;
            obj.sys = sys;
            
            obj.phyUBU = [0.3665, 6.1]; %% upper and lower bound limit
            obj.phyLBU = [-0.3665, 0];
            
            obj.defaultState = [12; -0.05; -0.05; -0.05]; % Valid for all state representations
            obj.defaultControl = [0];
            
            if nargin == 2
                obj.load_params_from_yaml(model_params_filepath);
            end
        end
        
        function [state_dot] = dynamics(obj, state, control, params)
            
            if nargin < 4
                params = [];
            end
            
            % Evaluate dynamics implementation for chosen state representation
            switch (obj.stateRep)
                
                case 'Flightpath'
                    state_dot = obj.dynamics_impl_flightpath(state, control, params);
                    
                case 'Pitch'
                    state_dot = obj.dynamics_impl_pitch(state, control, params);
                    
                case 'UW'
                    state_dot = obj.dynamics_impl_uw(state, control, params);
                    
                otherwise
                    error(['State representation ''' obj.stateRep ''' not found.']);
            end
            
        end
        
    end
    
    % Private methods
    methods (Access = protected)
        
        % Aerodynamic forces and moments implementation
        function [LIFT, DRAG, M] = lon_aero(obj, Va, alpha, q, dE)
            
            V0 = 12;
            
            CL = obj.CL0 + obj.CLa * alpha + obj.CLq * obj.c / (2.0 * V0) * q + obj.CLde * dE;
            CD = obj.CD0 + CL * CL / (pi * obj.e_oswald * obj.AR);
            Cm = obj.Cm0 + obj.Cma * alpha + obj.c / (2.0 * V0) * obj.Cmq * q + obj.Cmde * dE;
            
            dyn_press = 0.5 * obj.rho * Va * Va;
            LIFT = dyn_press * obj.S * CL;
            DRAG = dyn_press * obj.S * CD;
            M = dyn_press * obj.S * obj.c * Cm;
            
        end
        
        % Dynamics implementation for different state representations
        function [state_dot] = dynamics_impl_flightpath(obj, state, control, params)
            
            % Decompose state
            Va = state(1);
            alpha = state(2);
            q = state(3);
            gamma = state(4);
            
            dE = control(1);
            
            % Dynamics model
            [LIFT, DRAG, M] = obj.lon_aero(Va, alpha, q, dE);
            
            thrust = 0;
            b_thrust_ang = 0;
            M_thrust = 0;
            
            Va_dot = -DRAG / obj.mass + cos(alpha + b_thrust_ang) * thrust / obj.mass - obj.g * sin(gamma);
            gamma_dot = (LIFT / obj.mass + sin(alpha + b_thrust_ang) * thrust / obj.mass - obj.g * cos(gamma)) / Va;
            alpha_dot = q - gamma_dot;
            q_dot = (M + M_thrust) / obj.Iyy;
            
            % Compose state_dot
            state_dot = [Va_dot; alpha_dot; q_dot; gamma_dot];
            
        end
        
        function [state_dot] = dynamics_impl_pitch(obj, state, control, params)
            
            % Decompose state
            Va = state(1);
            alpha = state(2);
            q = state(3);
            pitch = state(4);
            
            dE = control(1);
            
            % Dynamics model
            [LIFT, DRAG, M] = obj.lon_aero(Va, alpha, q, dE);
            
            thrust = 0;
            b_thrust_ang = 0;
            M_thrust = 0;
            
            Va_dot = -DRAG / obj.mass + cos(alpha + b_thrust_ang) * thrust / obj.mass - obj.g * sin(pitch - alpha);
            alpha_dot = (-LIFT / obj.mass - sin(alpha + b_thrust_ang) * thrust / obj.mass + obj.g * cos(pitch - alpha) + q * Va) / Va;
            q_dot = (M + M_thrust) / obj.Iyy;
            pitch_dot = q;
            
            % Compose state_dot
            state_dot = [Va_dot; alpha_dot; q_dot; pitch_dot];
            
        end
        
        function [state_dot] = dynamics_impl_uw(obj, state, control, params)
            
            % Decompose state
            vx = state(1);
            vz = state(2);
            q = state(3);
            theta = state(4);
            
            dE = control(1);
            
            % Dynamics model
            Va = norm([vx, vz]);
            alpha = atan(vz / vx);
            [LIFT, DRAG, M] = obj.lon_aero(Va, alpha, q, dE);
            
            b_F_thrust = [0 0 0];
            M_thrust = 0;
            
            X = -cos(alpha) * DRAG + sin(alpha) * LIFT;
            Z = -cos(alpha) * LIFT - sin(alpha) * DRAG;
            vx_dot = (X + b_F_thrust(1)) / obj.mass - obj.g * sin(theta) - q * vz;
            vz_dot = (Z + b_F_thrust(3)) / obj.mass + obj.g * cos(theta) + q * vx;
            q_dot = (M + M_thrust) / obj.Iyy;
            theta_dot = q;
            
            % Compose state_dot
            state_dot = [vx_dot; vz_dot; q_dot; theta_dot];
            
        end
        
    end
    
end
