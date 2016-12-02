function [depth_cam,rgb_cam,Rdtrgb,Tdtrgb,image_name]=init(script)
    load(script);
    depth_cam=struct('K',K_ir,'DistCoef',distCoeffs_l);
    rgb_cam=struct('K',K_rgb,'DistCoef',distCoeffs_r);
    
    Rdtrgb=R12;
    Tdtrgb=T12;
    
    depth=dir('depth*.mat');
    images=dir('rgb*.png');
    for i=1:length(depth)    
        image_name(i)=struct('depth',depth(i).name,'rgb',images(i).name);
    end
 end