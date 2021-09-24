classdef TranscriptionBase
    %TRANSCRIPTIONBASE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        ocp
        opts
    end
    
    methods (Access = public)
        function [ocp, obj] = transcribe(obj, ocp, opts)
            obj.ocp = ocp;
            obj.opts = opts;
            
            % Create OptiStack using state and input dimenstions from OCP
            [ocp, obj] = obj.generate_decision_vars();
            obj.ocp = ocp;
            
            % Set cost using cost_impl from OCP
            obj.setup_cost();
            
            % Set dynamics constraints using cont_dynamics_impl from OCP
            obj.setup_dyn_constraints();
            
            % Set initial state
            ocp.opti.subject_to( obj.ocp.X_sym(:,1) == obj.ocp.x0_sym );
            
            % Return ocp
            obj.ocp = ocp;
            ocp = obj.ocp;
        end
        
        % Output functions
        function [U, T] = get_control_trajectory_impl(obj, Uopt, Topt, dt)
            [U, T] = obj.wrap_retime(Uopt, Topt, dt);
        end
        function [X, T] = get_state_trajectory_impl(obj, Xopt, Topt, dt)
            [X, T] = obj.wrap_retime(Xopt, Topt, dt);
        end
    end
    
    methods(Abstract)
        [ocp, obj] = generate_decision_vars(obj);
        setup_cost(obj);
        setup_dyn_constraints(obj);
    end
    
    methods (Static, Access = protected)
        function [X, T] = wrap_retime(X, Tx, dt)
            Xtt = timetable(seconds(Tx)', X');
            Xrt = retime(Xtt, 'regular', 'previous', 'TimeStep', seconds(dt));
            X = Xrt.Variables';
            T = seconds(Xrt.Time)';
        end
    end
end
