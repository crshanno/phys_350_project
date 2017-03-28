classdef BlenderMovePayload
    properties
        x;
        y;
        z;
    end
    methods
        function obj = BlenderMovePayload(x, y, z)
            obj.x = x;
            obj.y = y;
            obj.z = z;
        end
    end
end