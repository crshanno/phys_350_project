classdef BlenderCollisionPayload
    properties
        x;
        y;
        z;
        KE;
    end
    methods
        function obj = BlenderCollisionPayload(x, y, z, KE)
            obj.x = x;
            obj.y = y;
            obj.z = z;
            obj.KE = KE;
        end
    end
end