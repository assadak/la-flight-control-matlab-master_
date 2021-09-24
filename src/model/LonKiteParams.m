classdef LonKiteParams < handle
    
    properties (Access = protected)
        
        g = 9.806;
        rho = 1.1589;
        
        b
        c
        AR
        S
        
        mass
        Iyy
        
        e_oswald
        CD0
        
        CL0
        CLa
        Cm0
        Cma
        
        CLq
        Cmq
        
        CLde
        Cmde
        
    end
    
    methods (Access = public)
        
        function obj = load_params_from_yaml(obj, yaml_filepath)
            
            yaml = ReadYaml(yaml_filepath);
            
            obj.b = yaml.geom.b;
            obj.c = yaml.geom.c;
            obj.AR = yaml.geom.AR;
            obj.S = yaml.geom.S;
            
            obj.mass = yaml.inertia.mass;
            obj.Iyy = yaml.inertia.Iyy;
            
            obj.e_oswald = yaml.aero.e_oswald;
            obj.CD0 = yaml.aero.CD0;
            
            obj.CL0 = yaml.aero_aoa.CL0;
            obj.CLa = yaml.aero_aoa.CLa;
            
            obj.Cm0 = yaml.aero_aoa.Cm0;
            obj.Cma = yaml.aero_aoa.Cma;

            obj.CLq = yaml.aero_rate_pitch.CLq;
            obj.Cmq = yaml.aero_rate_pitch.Cmq;

            obj.CLde = yaml.aero_ctrl_elev.CLde;
            obj.Cmde = yaml.aero_ctrl_elev.Cmde;
            
        end

    end
    
end
