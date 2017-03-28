classdef BlenderCreatePayload
    properties
        x;
        y;
        z;
        radius;
    end
    methods
        function obj = BlenderCreatePayload(x, y, z, r)
            obj.x = x;
            obj.y = y;
            obj.z = z;
            obj.radius = r;
        end
    end
end