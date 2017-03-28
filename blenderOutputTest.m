function blenderOutputTest()    
    FILENAME = 'test.json';
    % initiate the output object
    b = BlenderOutput(FILENAME);
    
    % force flushing of buffers
    % (so that the file will be updated immediately after events are
    % triggered)
    b.setForceFlush(true);
    
    % create several objects
    for i = 1:50
        t = 0;
        x = 0;
        y = 0;
        z = 0;
        r = 1;
        b.objectCreate(i, t, x, y, z, r);
    end
    
    % move the objects over 10 seconds
    for t = 1:10
        for i = 1:50
            % all objects will move at different speeds
            x = t*i;
            y = t*i; 
            z = t*i;
            b.objectMove(i, t, x, y, z);
        end
    end
    
    % delete the objects
    for i = 1:50
        t = 10;
        b.objectDelete(i, t);
    end
    
    % close the output
    b.close();
end