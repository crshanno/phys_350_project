% Used to output data to be read by the blender visualisation
classdef BlenderOutput
   properties (Access = private)
       fileH;
       filename;
       forceFlush = false;
   end
   
   properties (Constant)
       defaultFilename = './log.json';
   end
   
   methods (Static)
       function clear(filename)
           % deletes the given file
           if nargin < 1
               filename = BlenderOutput.defaultFilename;
           end
          delete(filename);
       end
       
       function closeAll()
           % closes all open files
           fclose all;
       end
   end
   
   methods (Access = private)      
       function printJSON(obj, json)
           fprintf(obj.fileH, '%s\n', json);
           if (obj.forceFlush)
                drawnow('update'); % attempt to flush the write buffer
           end
       end 
   end
   
   methods (Access = public)
       function obj = BlenderOutput(filename)
           % initiates the output and returns an object for writing events
           % in the future
           if nargin < 1
               filename = BlenderOutput.defaultFilename;
           end
           obj.filename = filename;
           obj.fileH = fopen(obj.filename, 'a');
       end
     
       function close(obj)
           % closes the output. call this when the script is done
           fclose(obj.fileH);
       end
       
       function setForceFlush(obj, force)
           % set whether writen data should be flushed imediately.
           % useful if using blender in watch mode simultaneously
           obj.forceFlush = force;
       end
       
       function objectCreate(obj, idx, time, x, y, z, r)
           % trigger a creation event. ie. an object spawned
           
           % idx: the index of the object created
           % time: time in which it was created (seconds)
           % x, y, z: initial position, cartesian coordinates (meters)
           % r: radius (assuming sphere) (meters)
           b = BlenderFormat('CREATE', idx, time, BlenderCreatePayload(x, y, z, r));
           obj.printJSON(b.getJSON());         
       end
       
       function objectMove(obj, idx, time, x, y, z)
           % trigger a move event for the specfied object
           % these do not have to be frequent, since blender will
           % interpolate between move events, and too many will make the
           % import slow
           
           % idx: the index of the object moved
           % time: time in which it is at new position (seconds)
           % x, y, z: new position, cartesian coordinates (meters)
           b = BlenderFormat('MOVE', idx, time, BlenderMovePayload(x, y, z));
           obj.printJSON(b.getJSON());
       end
       
       function objectDelete(obj, idx, time)
           % trigger a delete event
           % the given object will be hidden from the render after the time
           % given
           
           % idx: the index of the object to be delete
           % time: time of deletion (seconds)
           b = BlenderFormat('DELETE', idx, time);
           obj.printJSON(b.getJSON());
       end
       
       function collision(obj, time, x, y, z, KE)
          % a collision event. 
          % currently unimplemented in the blender script, but it could be
          % used for audio/visual to emphasize the collision
          
          % time: time in which collision happened (seconds)
          % x, y, z: position of collision, cartesian coordinates (meters)
          % KE: kinetic energy of collision, used to scale effect (joules)
          b = BlenderFormat('COLLISION', 0, time, BlenderCollisionPayload(x, y, z, KE));
          obj.printJSON(b.getJSON());
       end
   end
end