classdef BlenderFormat
    properties (Access = public)
        event;
        idx;
        time;
        payload;
    end
    
    methods (Access = public)
        function obj = BlenderFormat(event, idx, time, payload)
            obj.event = event;
            obj.idx = idx;
            obj.time = time;
            if (nargin == 4)
                obj.payload = payload;
            end
        end
        
        function json = getJSON(obj)
            json = jsonencode(obj);
        end            
    end
end