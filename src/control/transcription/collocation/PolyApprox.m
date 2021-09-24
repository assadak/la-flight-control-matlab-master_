classdef PolyApprox
    properties (Access = public)
        tf          % Final time of trajectory
        Nseg        % Number of segments
        D           % degree of polynomial basis function
        
        Tsegs       % Start times of segments
        d_seg       % Duration of one segment = tf/Nseg
        Tau         % Nondimensional time points within a segment [0, 1]
    end
    
    methods (Access = public)
        function obj = PolyApprox(tf, Nseg, D)
            obj = obj.setup(tf, Nseg, D);
        end
        
        function approx_test(obj)
            
            Tseg = obj.Tsegs;
            
            % True signal: Sine scaled from 0 to Tf
            f_true = @(t) ( repmat( sin(t * 2*pi/obj.tf) , 2, 1) );
            df_true = @(t) ( repmat( 2*pi/obj.tf * cos(t * 2*pi/obj.tf), 2, 1) );
            
            % Set values at collocation points to true function values
            iSeg = 1;
            Xseg0 = f_true(Tseg(iSeg) + obj.d_seg * obj.Tau);
            
            % Evaluate from polynomial approximation
            xplot = linspace(Tseg(iSeg), Tseg(iSeg+1), 21);
            for i=1:length(xplot)
                tau_eval = (xplot(i) - Tseg(iSeg)) / obj.d_seg;
                yplot(:,i) = obj.approx(Xseg0, tau_eval);
                dyplot(:,i) = obj.dapprox(Xseg0, tau_eval);
            end
            
            close all
            figure
            x = linspace(0, obj.tf, 100);
            plot(x, f_true(x), x, df_true(x)); hold on
            plot(xplot, yplot, 'x', xplot, dyplot, 'x'); grid on
        end
    end
    
    methods (Access = public)
        function obj = setup(obj, tf, Nseg, D)
            obj.tf = tf;
            obj.Nseg = Nseg;
            obj.D = D;
            obj.d_seg = obj.tf/obj.Nseg;
            obj.Tsegs = 0:obj.d_seg:tf-obj.d_seg;
            obj.Tau = obj.cgl_points(D, 0, 1);
        end
        
        % Grid for polynomials
        function TauSeg = cgl_points(obj, D, a, b)
            D = D+1;
            js = 0:D;
            TauSeg = (b-a)/2 + (b-a)/2 * -cos(pi*js/D);
            TauSeg = TauSeg(1:end-1);
        end
        
        % Interpolate on segment using lagrange polynomials
        function x_i = approx(obj, ...
                X_seg, ... % State values at collocation points in this segment
                tau_eval)  % Absolute time of evaluation
            
            if size(X_seg,2) ~= obj.D+1
                error('Mismatch of provided values and poly degree');
            end
            
            x_i = 0;
            for jD = 0:obj.D % j = 0:obj.D
                x_i = x_i + obj.L(jD, tau_eval) * X_seg(:,jD+1);
            end
        end
        
        function dx_i = dapprox(obj, ...
                X_seg, ... % State values at collocation points in this segment
                tau_eval)  % Absolute time of evaluation
            
            if size(X_seg,2) ~= obj.D+1
                error('Mismatch of provided values and poly degree');
            end
            
            dx_i = 0;
            for jD = 0:obj.D
                dx_i = dx_i + obj.dL(jD, tau_eval) * X_seg(:,jD+1);
            end
            dx_i = 1/obj.d_seg * dx_i;
        end
    end
    
    methods (Access = public)
        % Lagrange polynomial
        function L = L(obj, j, tau_eval)
            tau = obj.Tau;
            L = 1;
            for kD = 0:obj.D
                if kD ~= j
                    L = L * (tau_eval - tau(kD+1)) / (tau(j+1) - tau(kD+1));
                end
            end
        end
        
        % First derivative of Lagrange polynomial
        function dL = dL(obj, j, tau_eval)
            tau = obj.Tau;
            
            dL = 0;
            for iD = 0:obj.D
                if iD ~= j
                    % This is the product
                    product = 1;
                    for kD = 0:obj.D
                        if (kD ~= iD && kD ~= j)
                            product = product * (tau_eval - tau(kD+1)) / (tau(j+1) - tau(kD+1));
                        end
                    end % ---------------
                    dL = dL + product / (tau(j+1) - tau(iD+1));
                end
            end
        end
        
    end
end

