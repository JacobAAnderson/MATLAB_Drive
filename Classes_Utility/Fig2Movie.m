% Record frames and make movie


classdef Fig2Movie
    
    properties
        Frames
        Name
        Idx
        FrameRate
        FrameSize
    end
    
    methods
        
        function obj = Fig2Movie(name, sz), warning on backtrace
            
            obj.Frames = struct('cdata',[],'colormap',[]);
            obj.Name = name;
            obj.Idx = 1;
            obj.FrameRate = 30;
            obj.FrameSize = [];
            
            if nargin > 1                                                   % Prealocate memory if size is provided
                obj.Frames(sz) = struct('cdata',[],'colormap',[]);
            end
            
        end
        
        
        function obj = Add_Frame(obj, h)
            
            obj.Frames(obj.Idx) = getframe(h);
            obj.Idx = obj.Idx + 1;
            
        end
        
        
        function obj = Add_StillFrame(obj, h, sec)
            
            num = sec * obj.FrameRate;
            
            for ii = 1:num
                obj.Frames(obj.Idx) = getframe(h);
                obj.Idx = obj.Idx + 1;
            end
            
        end
        
        
        function MakeVideo(obj, path, frameRate)
            
            sz = cellfun(@size, {obj.Frames.cdata}, 'UniformOutput', false);
            sz = cat(1, sz{:});
            
            M = median(sz,1);
            
            fprintf("\n\nWritting Video to:\n%s\n\n",path)
            
            v = VideoWriter(path);
            v.Quality = 100;
            
            if nargin > 2, v.FrameRate = frameRate; end
            
            open(v)
            
            try
                for ii = 1: size(obj.Frames,2)
                    
                    f = obj.Frames(ii);
                    
                    if any(size(f.cdata) ~= M) 
                        f.cdata = imresize(f.cdata, [M(1), M(2)]); 
                    end 
                    
                    writeVideo(v, f);
                end
                
                close(v)
                
            catch ME
                close(v)
                rethrow(ME)
                
            end
            
        end
        
        
    end
end




