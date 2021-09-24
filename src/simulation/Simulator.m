classdef Simulator < handle
    
    properties
        dt
        dyn_func
    end
    
    methods (Access = public)
        
        function obj = Simulator(dt, dyn_func)
            
            obj.dt = dt;
            
            if nargin == 2
                obj.dyn_func = dyn_func;
            end
        end
        
        function xp = simulate_step(obj, x0, u, p, dyn_func, dt)
            
            if nargin < 4
                p = [];
            end
            if nargin < 5
                dyn_func = obj.dyn_func;
            end
            if nargin < 6
                dt = obj.dt;
            end
            
            xp = obj.simulate_nonlinear_step(dyn_func, x0, u, p, dt);
        end
        
        function X = simulate(obj, x0, U, p, dyn_func, dt)
            
            if nargin < 4
                p = [];
            end
            if nargin < 5
                dyn_func = obj.dyn_func;
            end
            if nargin < 6
                dt = obj.dt;
            end
            
            X = obj.simulate_nonlinear(dyn_func, dt, x0, U, p);
        end
        
        function plot_state_trajectory(obj, time, state_trajs, sys)

            % Prepare cell array of structs
            if iscell(state_trajs)
                if iscell(state_trajs(1))
                    % Two-level cell, most general case, do nothing
                else
                    % One-level cell
                    state_trajs = {state_trajs};
                end
            else
                % No cell
                state_trajs = {{state_trajs, 'traj'}};
            end
            
            % Work with generalized cell array
            nTraj = length(state_trajs);
            nx = size(state_trajs{1}{1}, 1);
            
            figure;
            convert_rad = false;
            
            for iTraj = 1:nTraj
                for ix = 1:nx
                    subplot(nx, 1, ix);
                    
                    if nargin > 3
                        convert_rad = contains(sys.StateUnit{ix}, 'rad');
                        if convert_rad
                            plot(time, state_trajs{iTraj}{1}(ix,:));
                            strrep(sys.StateUnit{ix}, 'rad', 'deg');
                        else
                            plot(time, state_trajs{iTraj}{1}(ix,:));
                            ylabel(sys.StateUnit{ix});
                        end
                        title(sys.StateName{ix});
                    else
                        plot(time, state_trajs{iTraj}{1}(ix,:));
                        
                    end
                    hold on; grid on;
                    
                end
                subplot(nx, 1, 1);
                
                legendEntries{iTraj} = state_trajs{iTraj}{2};
            end
            
            legend(legendEntries);
        end
        
        function ph = plot_control_state_trajectory(obj, time, control_traj, state_trajs, sys)
            % Prepare cell array of structs
            if iscell(state_trajs)
                if iscell(state_trajs(1))
                    % Two-level cell, most general case, do nothing
                else
                    % One-level cell
                    state_trajs = {state_trajs};
                end
            else
                % No cell
                state_trajs = {{state_trajs, 'traj'}};
            end
            
            % Work with generalized cell array
            nTraj = length(state_trajs);
            nx = size(state_trajs{1}{1}, 1);
            nu = size(control_traj, 1);
            nxu = nx + nu;
            
            figure;
            convert_rad = false;
            
            for iu = 1:nu
                ph.ax(iu) = subplot(nxu, 1, iu);
                
                if nargin > 4
                    convert_rad = contains(sys.InputUnit{iu}, 'rad');
                    if convert_rad
                        plot(time(1:end-1), 180/pi * control_traj(iu,:));
                        ylabel( strrep(sys.InputUnit{iu}, 'rad', 'deg') );
                    else
                        plot(time(1:end-1), control_traj(iu,:));
                        ylabel(sys.InputUnit{iu});
                    end
                    title(sys.InputName{iu});
                else
                    plot(time(1:end-1), control_traj(iu,:));
                end
                hold on; grid on;
                YLIM = ylim();
                ylim(1.1*YLIM)
                
            end
            
            for iTraj = 1:nTraj
                for ix = 1:nx
                    ph.ax(nu+ix) = subplot(nxu, 1, nu+ix);
                    
                    if nargin > 4
                        % System properties are available
                        convert_rad = contains(sys.StateUnit{ix}, 'rad');
                        
                        if convert_rad
                            plot(time, 180/pi * state_trajs{iTraj}{1}(ix,:));
                            ylabel( strrep(sys.StateUnit{ix}, 'rad', 'deg') );
                        else
                            plot(time, state_trajs{iTraj}{1}(ix,:));
                            ylabel( sys.StateUnit{ix} );
                        end
                        title(sys.StateName{ix});
                    else
                        plot(time, state_trajs{iTraj}{1}(ix,:));
                    end
                    hold on; grid on;
                    
                end
                subplot(nxu, 1, nu+1);
                
                legendEntries{iTraj} = state_trajs{iTraj}{2};
            end
            
            legend(legendEntries);
        end
        
    end
    
    methods (Access = protected)
        
        % Stepper for nonlinear dynamics
        function [x] = simulate_nonlinear_step(obj, dyn_func, x0, u, p, dt)
            % Integrate nonlinear system for dt
            test=1;
            [~, X_] = ode45( @(t_, x_) dyn_func(x_, u, p), [0 dt], x0);
            x = X_(end,:)';
        end
        
        % Integration along input trajectory
        function [X] = simulate_nonlinear(obj, dyn_func, dt, x0, U, p)
            
            nSteps = size(U, 2);
            X = [x0, zeros(size(x0, 1), nSteps)];
            
            for iStep = 1:nSteps
                X(:, iStep+1)  = obj.simulate_nonlinear_step(dyn_func, X(:,iStep), U(:,iStep), p, dt);
            end
        end
        
    end
    
end
