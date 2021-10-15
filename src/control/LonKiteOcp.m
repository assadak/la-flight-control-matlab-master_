classdef LonKiteOcp < OcpBase
    
    methods (Static, Access = public)
        function control_objectives = ControlObjectives()
            control_objectives = struct(...
                'track_angle',      1, ...
                'track_airspeed',   0, ...
                'minimize_control', 1, ...
                'track_height',    0);
        end
    end
    
    properties (Access = protected)
        control_objectives
    end
    
    methods (Access = public)    %% dynamics, objective,end time, cont options, ?
        function obj = LonKiteOcp(lonKite, control_objectives, tf, opts, varargin)
            
            obj@OcpBase(lonKite, tf);
            
            obj.control_objectives = control_objectives;
            
            % Parse additional function parameters
            p = inputParser;
            p.addOptional('transcription', 'multiple_shooting');
            p.addOptional('use_cos_grid', 0, @islogical);
            p.addOptional('N', opts.N, @isnumeric);
            p.addOptional('approx_control', 1, @islogical);
            p.parse(varargin{:});
            
            transcript_method = p.Results.transcription;
            opts.use_cos_grid = p.Results.use_cos_grid;
            opts.N = p.Results.N;
            opts.approx_control = p.Results.approx_control;
            
            obj.opts = opts;
            
            % Create symbolic parameters (used in cost implementation)
            obj.p_sym.config.angle_ref    = obj.opti.parameter(1, 1);  % Reference angle
            obj.p_sym.config.Va_ref       = obj.opti.parameter(1, 1);  % Reference airspeed
            obj.p_sym.config.h_ref        = obj.opti.parameter(1, 1);  % Reference height
            
            obj.p_sym.tuning.mayer_multiplier = obj.opti.parameter(1, 1);  % Multiplier for non_control stage cost to obtain mayer cost
            
            obj.p_sym.tuning.W_angle_err  = obj.opti.parameter(1, 1);  % Weight on angle error
            obj.p_sym.tuning.W_Va_err     = obj.opti.parameter(1, 1);  % Weight on airspeed error
            obj.p_sym.tuning.W_h_err      = obj.opti.parameter(1, 1);  % Weight on height error
            
            obj.p_sym.tuning.R_diag       = obj.opti.parameter(lonKite.nu, 1);  % Diagonal of control weight matrix
            
            % Setup optimization problem
            obj = obj.setup(transcript_method);
        end
        
        % Override functions from base class ------------------------------
        function cost = stage_cost_impl(obj, x, u)
            cost = obj.get_non_control_cost(x) + obj.get_control_cost(u);
        end
        
        function cost = terminal_cost_impl(obj, x)
            cost = obj.p_sym.tuning.mayer_multiplier * obj.get_non_control_cost(x);
        end
        
        function x_dot = cont_dynamics_impl(obj, x, u)
            x_dot = obj.lonKite.dynamics(x, u);
        end
        
        % OCP-specific setters --------------------------------------------
        % Setters for config parameters
        function set_angle_ref(obj, angle_ref)
            obj.opti.set_value( obj.p_sym.config.angle_ref, angle_ref);
        end
        function set_Va_ref(obj, Va_ref)
            obj.opti.set_value( obj.p_sym.config.Va_ref, Va_ref);
        end
        function set_h_ref(obj, h_ref)
            obj.opti.set_value( obj.p_sym.config.h_ref, h_ref);
        end
        
        % Setters for tuning parameters
        function set_mayer_multiplier(obj, mayer_multiplier)
            obj.opti.set_value( obj.p_sym.tuning.mayer_multiplier, mayer_multiplier);
        end
        
        function set_W_Va_err(obj, W_Va_err)
            obj.opti.set_value( obj.p_sym.tuning.W_Va_err, W_Va_err);
        end
        function set_W_angle_err(obj, W_angle_err)
            obj.opti.set_value( obj.p_sym.tuning.W_angle_err, W_angle_err);
        end
        function set_W_h_err(obj, W_h_err)
            obj.opti.set_value( obj.p_sym.tuning.W_h_err, W_h_err);
        end
        
        function set_R_diag(obj, R_diag)
            obj.opti.set_value( obj.p_sym.tuning.R_diag, R_diag);
        end
    end
    
    methods (Access = private)
        % Helpers for non-/control stage costs
        function non_control_cost = get_non_control_cost(obj, x)
            
            non_control_cost = 0;
            
            if obj.control_objectives.track_angle
                % Objective: Track angle
                angle_err = obj.p_sym.config.angle_ref - x(4);
                non_control_cost = non_control_cost + obj.p_sym.tuning.W_angle_err * angle_err^2;
            end
            
            if obj.control_objectives.track_airspeed
                % Objective: Track airspeed
                Va_err = obj.p_sym.config.Va_ref - x(1);
                non_control_cost = non_control_cost + obj.p_sym.tuning.W_Va_err * Va_err^2;
            end
            
            if obj.control_objectives.track_height
                height_err = obj.p_sym.config.h_ref -(-x(6));
                non_control_cost = non_control_cost + obj.p_sym.tuning.W_h_err * height_err^2;
            end
            
            
        %% add my code here?
% does it itirate over this file or it just takes the difference of the
% whole vector at once?

            
            
            
            
            
            
            
        end
        
        function control_cost = get_control_cost(obj, u)
            
            control_cost = 0;
            
            if obj.control_objectives.minimize_control
                % Objective: Minimize control effort
                control_cost = control_cost + u' * diag(obj.p_sym.tuning.R_diag) * u;
            end
        end
    end
    
end

