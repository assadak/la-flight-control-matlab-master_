classdef KiteDynamics < KiteParams
    
    properties (Access = public)
        
        nx = 13;         % Number of states
        nu = 3;         % Number of inputs
        
        defaultState    % State for which the dynamics is well defined
        defaultControl  % Control for which the dynamics is well defined
        phyUBU          % Physical upper control bound
        phyLBU          % Physical Lower control bound
        
        sys             % Struct containing system properties
        
    end
    
    % Public methods
    methods (Access = public)
        
        function obj = KiteDynamics(model_params_filepath)
            
            sys.StateName = {'vx' 'vy' 'vz' 'wx' 'wy' 'wz' 'posN' 'posE' 'posD' 'qw' 'qx' 'qy' 'yz'};
            sys.StateUnit = {'m/s' 'm/s' 'm/s' 'rad/s' 'rad/s' 'rad/s' 'm' 'm' 'm' '-' '-' '-' '-'};
            sys.InputName = {'dE' 'dR' 'dA'};
            sys.InputUnit = {'rad' 'rad' 'rad'};
            sys.OutputName = sys.StateName;
            sys.OutputUnit = sys.StateUnit;
            obj.sys = sys;
            
            obj.phyUBU = [0.3665, 0.4625, 0.3578];
            obj.phyLBU = [-0.3665, -0.4625, -0.3578];
            
            obj.defaultState = [12; 0; 0; ...
                0; 0; 0; ...
                0; 0; -100; ...
                1; 0; 0; 0];
            obj.defaultControl = [0; 0; 0];
            
            obj.load_params_from_yaml(model_params_filepath);
        end
        
        function [state_dot] = dynamics(obj, state, control, params)
            
            if nargin < 4
                params = [];
            end
            
            state_dot = obj.dynamics_impl(vec(state), vec(control), vec(params));
            
        end
        
        function ph = plot_pose_trajectory(obj, X, name, ph)
            
            if nargin < 4
                ph = [];
            end
            
            if nargin < 3
                name = 'traj';
            end
            
            if size(X, 2) > size(X, 1)
                X = X';
            end
            dat = array2table(X, 'VariableNames', obj.sys.StateName);
            
            figTitle = 'Local position';
            if isempty(ph)
                f = figure;
                f.Name = figTitle;
                set(0, 'currentfigure', f);
                clf;
                ax = subplot(1, 1, 1);
            else
                f = ph.fig;
                ax = ph.ax;
            end
            
            set(f, 'currentaxes', ax(1));
            % Position trajectory
            p = plot3(dat.posE, dat.posN, dat.posD, 'b'); hold on
            if ~isempty(name), p.DisplayName = name; end
            hold on;
            
            % Discretized attitude
            p = plot3(dat.posE(1), dat.posN(1), dat.posD(1), '^k');
            p.Annotation.LegendInformation.IconDisplayStyle = 'off';
            legend('-DynamicLegend');
            
            ax(1).ZDir = 'reverse';
            
            title(figTitle)
            ylabel('North')
            xlabel('East')
            zlabel('Down')
            
            grid on
            
            axis equal
            
            yLim = ylim;
            xLim = xlim;
            zLim = zlim;
            
            ySpan = diff(yLim);
            xSpan = diff(xLim);
            zSpan = diff(zLim);
            
            ylim([yLim(1) - 0.1 * max(ySpan, 20), yLim(2) + 0.1 * max(ySpan, 20)]);
            xlim([xLim(1) - 0.1 * max(xSpan, 20), xLim(2) + 0.1 * max(xSpan, 20)]);
            zlim([zLim(1) - 0.1 * max(zSpan, 20), zLim(2) + 0.1 * max(zSpan, 20)]);
            
            ph.fig = f;
            ph.ax = ax;
        end
        
    end
    
    % Private methods
    methods (Access = protected)
        function [state_dot] = dynamics_impl(obj, state, control, params)
            
            %% TODO:
            %Insert 'obj.' before each parameter (defined in KiteParams.m),
            %Translate tether model.
            %%
            
            % Decompose state
            v = state(1:3);
            w = state(4:6);
            pos = state(7:9);
            q = state(10:13)';
            q_bn = quatinv(q);
            
            dE = control(1);
            dR = control(2);
            dA = control(3);
            
            n_vW = params(1:3);
            
            % Dynamics model
            %%% Aerodynamic (Wind) frame %%%
            % Wind and apparent velocity in body frame
            b_vW = obj.quatTransform(q_bn, n_vW);
            b_va = v - b_vW;
            
            % Aerodynamic variables (Airspeed, angle of attack, side slip angle)
            Va = norm(b_va);
            alpha = atan(b_va(3) / b_va(1));
            beta = asin(b_va(2) / Va);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% Aerodynamic Forces and Moments in aerodynamic (wind) frame %%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            q_ba = quatmultiply(obj.T2quat(alpha), obj.T3quat(-beta));
            
            % TODO: Implement option for  V0 = (FIX_V0) ? FIXED_V0 : Va;
            V0 = Va;
            
            dyn_press = 0.5 * obj.rho * Va^2;
            
            CL = obj.CL0 + obj.CLa * alpha + obj.CLq * obj.c / (2.0 * V0) * w(1) + obj.CLde * dE;
            CD = obj.CD0 + CL^2 / (pi * obj.e_oswald * obj.AR);
            
            % Forces in x, y, z directions: -Drag, Side force, -Lift
            LIFT = dyn_press * obj.S * CL;
            DRAG = dyn_press * obj.S * CD;
            SF = dyn_press * obj.S * (obj.CYb * beta + obj.b / (2.0 * V0) * (obj.CYp * w(1) + obj.CYr * w(3)) + obj.CYdr * dR);
            
            a_F_aero = [-DRAG; SF; -LIFT];
            
            % Moments about x, y, z axes: L, M, N
            L = dyn_press * obj.S * obj.b * (obj.Cl0 + obj.Clb * beta + obj.b / (2.0 * V0) * (obj.Clp * w(1) + obj.Clr * w(3)) + obj.Clda * dA + obj.Cldr * dR);
            
            M = dyn_press * obj.S * obj.c * (obj.Cm0 + obj.Cma * alpha + obj.c / (2.0 * V0) * obj.Cmq * w(2) + obj.Cmde * dE);
            
            N = dyn_press * obj.S * obj.b * (obj.Cn0 + obj.Cnb * beta + obj.b / (2.0 * V0) * (obj.Cnp * w(1) + obj.Cnr * w(3)) + obj.Cnda * dA + obj.Cndr * dR);
            
            b_M_aero = [L; M; N];
            
            % Aerodynamic Forces and Moments in body frame
            b_F_aero = obj.quatTransform(q_ba, a_F_aero);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% Graviation, Tether (body frame %%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Gravitational acceleration
            b_g = obj.quatTransform(q_bn, [0; 0; obj.g]);
            
            % Tether force and moment
            b_F_tether = [0;0;0];
            %     if (teth_on)
            %     {
            %         SX dist = SX::norm_2(pos);
            %         const double tethLen = 120.0;
            %         const double tethDiameter = 0.3e-3;
            %         const double tethCrossArea = M_PI_4 * tethDiameter * tethDiameter;
            %         const double tethDensity = 1.15;
            %         const double tethMassPerMeter = tethCrossArea * 1.0 * tethDensity;
            %         const double tethE = 9.05e9;
            %
            %         const double c_orth = 1.2;
            %
            %         /* Weight_tether */
            %         //SX decl = SX::atan(SX::norm_2(r(Slice(0, 2))) / -r(2));
            %         //SX g_Wteth = 0.5 * (1 + cos(decl)) * tethLen * tethMassPerMeter * SX::vertcat({0, 0, g}); // probably invalid
            %         SX n_Wteth = tethLen * tethMassPerMeter * SX::vertcat({0, 0, g});
            %         SX b_Wteth = polymath::quat_transform(q_bn, n_Wteth);
            %
            %         /* Drag_tether */
            %         SX q_ba = polymath::quat_multiply(polymath::T2quat(alpha), polymath::T3quat(-beta));
            %
            %         SX lat = SX::atan2(-pos(0), SX::norm_2(pos(Slice(1, 3))));
            %         SX lon = SX::atan2(-pos(1), -pos(2));
            %         SX q_ln = polymath::quat_multiply(polymath::T2quat(-lat), polymath::T1quat(lon));
            %         SX q_nl = polymath::quat_inverse(q_ln);
            %         SX q_lb = polymath::quat_multiply(q_ln, q);
            %
            %         SX l_va = polymath::quat_transform(q_lb, b_va);
            %         SX l_va_proj = SX::vertcat({l_va(0), l_va(1), 0});
            %
            %         SX q_bl = polymath::quat_multiply(q_bn, q_nl);
            %         SX f_va_proj = polymath::quat_transform(q_bl, l_va_proj);
            %
            %         SX b_Dteth = 0.125 * tethDensity * -f_va_proj * SX::norm_2(f_va_proj) * c_orth * tethDiameter * tethLen;
            %
            %         /* LonForce_tether */
            %         SX tethElongation = (dist - tethLen) / tethLen;
            % //        SX g_lonFteth = -pos / dist * SX::fmax(0, (tethE * tethCrossArea) * (tethElongation + 0.005));
            %         SX g_lonFteth = -pos / dist * 5;
            %         SX b_lonFteth = polymath::quat_transform(q_bn, g_lonFteth);
            %         b_F_tether = b_lonFteth + b_Dteth + b_Wteth;
            %     }
            %     else
            %         b_F_tether = SX::vertcat({0, 0, 0});
            %
            tether_outlet_position = [obj.teth_rx; obj.teth_ry; obj.teth_rz];
            b_M_tether = cross(tether_outlet_position, b_F_tether);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% Motion equations (body frame) %%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Linear motion equation
            spec_nongrav_force = (b_F_aero + b_F_tether) / obj.mass;
            v_dot = spec_nongrav_force + b_g - cross(w, v);
            
            % Angular motion equation
            J = diag([obj.Ixx, obj.Iyy, obj.Izz]);
            w_dot = inv(J) * (b_M_aero + b_M_tether - cross(w, J*w));
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% Kinematic Equations (geodetic frame) %%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Translation: Aircraft position derivative
            r_dot = obj.quatTransform(q, v);                            % q = q_nb, v (body frame)
            
            % Rotation: Aircraft attitude derivative
            % Quaternion representation
            lambda = -5;
            q_dot = 0.5 * quatmultiply(q, [0, w']) ...              % q = q_nb, w = omega (body frame)
                + 0.5 * lambda * q * (dot(q, q) - 1);        % Quaternion norm stabilization term,
            % as in Gros: 'Baumgarte Stabilisation over the SO(3) Rotation Group for Control',
            % improved: lambda negative and SX::dot(q, q) instead of lambda positive and 1/SX::dot(q, q).
            
            % Compose state_dot
            state_dot = [v_dot; w_dot; r_dot; q_dot'];
        end
    end
    
    methods (Static, Access = private)
        function q = T1quat(rotAng)
            q = [cos(-rotAng / 2.0), sin(-rotAng / 2.0), 0, 0];
        end
        function q = T2quat(rotAng)
            q = [cos(-rotAng / 2.0), 0, sin(-rotAng / 2.0), 0];
        end
        function q = T3quat(rotAng)
            q = [cos(-rotAng / 2.0), 0, 0, sin(-rotAng / 2.0)];
        end
        
        function b_vec = quatTransform(q_ba, a_vec)
            if size(a_vec,1) == 3
                % Quat is 1x4, Vector is 3x1
                res_temp = quatmultiply( q_ba, quatmultiply([1; a_vec]', quatinv(q_ba)) );
                b_vec = res_temp(:, 2:4)';
                
            elseif size(a_vec,2) == 3
                % Quat is nx4, Vector is nx3
                res_temp = quatmultiply( q_ba, quatmultiply([ones(size(a_vec,1),1), a_vec], quatinv(q_ba)) );
                b_vec = res_temp(:, 2:4);
            end
            
        end
    end
    
end
