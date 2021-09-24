classdef KiteParams < handle
    
    properties (Access = protected)
        
        g = 9.806;
        rho = 1.1589;
        
        %% Airplane geometry
        b  % Wing span
        c  % Chord length
        AR % Aspect ratio
        S  % Wing surface
        
        %% Mass and inertia
        mass
        Ixx
        Iyy
        Izz
        
        %% Aerodymic effects
        % Drag
        e_oswald % Oswald efficiency
        CD0      % Minimum drag (at CL = 0)
        
        % ... from angle of attack
        CL0
        CLa % Lift ...
        Cm0
        Cma % Pitch moment ...
        
        % ... from side slip
        CYb % Side force ...
        Cl0
        Clb % Roll moment ...
        Cn0
        Cnb % Yaw moment ...
        
        % ... from pitch rate
        CLq % Lift ...
        Cmq % Pitch moment ...
        
        % ... from roll rate
        CYp % Side force ...
        Clp % Roll moment ...
        Cnp % Yaw moment ...
        
        % ... from yaw rate
        CYr % Side force ...
        Clr % Roll moment ...
        Cnr % Yaw moment ...
        
        %% Aerodynamic effects of controls
        % ... from elevator
        CLde % Lift ...
        Cmde % Pitch moment ...
        
        % ... from ailerons
        Clda % Roll rate ...
        Cnda % Yaw rate ...
        
        % ... from rudder
        CYdr % Side force ...
        Cldr % Roll moment ...
        Cndr % Yaw moment ...
        
        %% Tether parameters
        % Attachment point in aircraft body frame
        teth_rx
        teth_ry
        teth_rz
        
        teth_length % Nominal tether length
        teth_Ks     % Spring coefficient
        teth_Kd     % Damping coefficient
        
    end
    
    methods (Access = public)
        
        function obj = load_params_from_yaml(obj, yaml_filepath)
            
            yaml = ReadYaml(yaml_filepath);
            
            obj.b = yaml.geom.b;
            obj.c = yaml.geom.c;
            obj.AR = yaml.geom.AR;
            obj.S = yaml.geom.S;
            
            obj.mass = yaml.inertia.mass;
            obj.Ixx = yaml.inertia.Ixx;
            obj.Iyy = yaml.inertia.Iyy;
            obj.Izz = yaml.inertia.Izz;
            
            obj.e_oswald = yaml.aero.e_oswald;
            obj.CD0 = yaml.aero.CD0;
            
            obj.CL0 = yaml.aero_aoa.CL0;
            obj.CLa = yaml.aero_aoa.CLa;
            obj.Cm0 = yaml.aero_aoa.Cm0;
            obj.Cma = yaml.aero_aoa.Cma;
            
            obj.CYb = yaml.aero_ss.CYb;
            obj.Cl0 = yaml.aero_ss.Cl0;
            obj.Clb = yaml.aero_ss.Clb;
            obj.Cn0 = yaml.aero_ss.Cn0;
            obj.Cnb = yaml.aero_ss.Cnb;
            
            obj.CLq = yaml.aero_rate_pitch.CLq;
            obj.Cmq = yaml.aero_rate_pitch.Cmq;
            
            obj.CYp = yaml.aero_rate_roll.CYp;
            obj.Clp = yaml.aero_rate_roll.Clp;
            obj.Cnp = yaml.aero_rate_roll.Cnp;
            
            obj.CYr = yaml.aero_rate_yaw.CYr;
            obj.Clr = yaml.aero_rate_yaw.Clr;
            obj.Cnr = yaml.aero_rate_yaw.Cnr;
            
            obj.CLde = yaml.aero_ctrl_elev.CLde;
            obj.Cmde = yaml.aero_ctrl_elev.Cmde;
            
            obj.Clda = yaml.aero_ctrl_ail.Clda;
            obj.Cnda = yaml.aero_ctrl_ail.Cnda;
            
            obj.CYdr = yaml.aero_ctrl_rud.CYdr;
            obj.Cldr = yaml.aero_ctrl_rud.Cldr;
            obj.Cndr = yaml.aero_ctrl_rud.Cndr;
            
            
            obj.teth_rx = yaml.tether.teth_rx;
            obj.teth_ry = yaml.tether.teth_ry;
            obj.teth_rz = yaml.tether.teth_rz;
            
            obj.teth_length = yaml.tether.teth_length;
            obj.teth_Ks = yaml.tether.teth_Ks;
            obj.teth_Kd = yaml.tether.teth_Kd;
        end
        
    end
    
end
