classdef LonKiteDynamicsWithPos < LonKiteDynamics
    
    % Public methods
    methods (Access = public)
        
        function obj = LonKiteDynamicsWithPos(state_rep, model_params_filepath)
            
            obj@LonKiteDynamics(state_rep, model_params_filepath);
            
            obj.nx = obj.nx + 2;
            
            sys = obj.sys;
            sys.StateName = [obj.sys.StateName, {'posHor', 'posDown'}];
            sys.StateUnit = [obj.sys.StateUnit, {'m', 'm'}];
            
            sys.OutputName = obj.sys.StateName;
            sys.OutputUnit = obj.sys.StateUnit;
            obj.sys = sys;
           
            obj.defaultState = [obj.defaultState; 0; -100]; % Valid for all state representations
            
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
    
    % Ptivate methods
    methods (Access = protected)
        
        % Dynamics implementation for different state representations
        function [state_dot] = dynamics_impl_flightpath(obj, state, control, params)
            
            % Decompose state
            core_state = state(1:4);
            core_control = control;
            
            Va = state(1);
            gamma = state(4);
            
            % Dynamics model
            core_state_dot = dynamics_impl_flightpath@LonKiteDynamics(obj, core_state, core_control, params);
            
            vHor = Va * cos(gamma);
            vDown = -Va * sin(gamma);
            
            % Compose state_dot
            state_dot = [core_state_dot; vHor; vDown];
            
        end
        
        function [state_dot] = dynamics_impl_pitch(obj, state, control, params)
            
            % Decompose state
            core_state = state(1:4);
            core_control = control;
            
            Va = state(1);
            gamma = state(4) - state(2); % pitch - alpha
            
            % Dynamics model
            core_state_dot = dynamics_impl_pitch@LonKiteDynamics(obj, core_state, core_control, params);
            
            vHor = Va * cos(gamma);
            vDown = -Va * sin(gamma);
            
            % Compose state_dot
            state_dot = [core_state_dot; vHor; vDown];
            
        end
        
        function [state_dot] = dynamics_impl_uw(obj, state, control, params)
            
            % Decompose state
            core_state = state(1:4);
            core_control = control;
            
            vx = state(1);
            vz = state(2);
            pitch = state(4);
            
            % Dynamics model
            core_state_dot = dynamics_impl_uw@LonKiteDynamics(obj, core_state, core_control, params);
            
            vHor = vx * cos(pitch) + vz * sin(pitch);
            vDown = -vx * sin(pitch) + vz * cos(pitch);
            
            % Compose state_dot
            state_dot = [core_state_dot; vHor; vDown];
            
        end
        
    end
    
end
