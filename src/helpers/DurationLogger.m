classdef DurationLogger
    
    properties (Access = public)
        storage
    end
    
    properties (Access = public)
        time_unit_in = 'Seconds';
        time_unit_out = 'MilliSeconds'
    end
    
    methods (Access = public)
        function obj = log(obj, value)
            obj.storage = [obj.storage, value];
        end
        
        function printStatistics(obj)
            if (obj.time_unit_out == "MicroSeconds")
                unit_out_str = 'us';
                if (obj.time_unit_in == "MicroSeconds")
                    scaling = 1;
                elseif (obj.time_unit_in == "MilliSeconds")
                    scaling = 1e3;
                else
                    scaling = 1e6;
                end
            elseif (obj.time_unit_out == "Seconds")
                unit_out_str = 's';
                if (obj.time_unit_in == "MicroSeconds")
                    scaling = 1e-6;
                elseif (obj.time_unit_in == "MilliSeconds")
                    scaling = 1e-3;
                else
                    scaling = 1;
                end
            else % time_unit_out == MilliSeconds
                unit_out_str = 'ms';
                if (obj.time_unit_in == "MicroSeconds")
                    scaling = 1e-3;
                elseif (obj.time_unit_in == "MilliSeconds")
                    scaling = 1;
                else % time_unit_in == Seconds
                    scaling = 1e3;
                end
            end
            
            max_val = max(obj.storage);
            min_val = min(obj.storage);
            median_val = median(obj.storage);
            
            out_str = ['\nComputation statistics [', unit_out_str, ']: ' ...
                '\nmax: ' num2str(max_val * scaling) ...
                '\nmed: ' num2str(median_val * scaling) ...
                '\nmin: ' num2str(min_val * scaling) ...
                '\n'];
            fprintf(out_str);
        end
        function obj = reset(obj)
            storage = [];
        end
    end
    
end

