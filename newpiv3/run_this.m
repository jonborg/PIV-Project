close all
clear all

script='camParams.mat';
[depth_cam,rgb_cam,Rdtrgb,Tdtrgb,image_name]=init(script);
[pcloud, transforms]=reconstruction(image_name, depth_cam, rgb_cam, Rdtrgb,Tdtrgb);
