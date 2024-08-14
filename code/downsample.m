%% Downsample to coarser resolution 
function [udown,tdown] = downsample(u,t,down_fps,frame_rate)
    delf = frame_rate/down_fps;
    if (floor(delf)~=delf) 
        disp(['Error delf = frame_rate/down_fps is not an integer -->' num2str(delf)]);
    end
    tdown = t(1:delf:end); 
    udown = u(1:delf:end);
end


