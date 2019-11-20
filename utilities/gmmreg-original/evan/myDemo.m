
% mydemo
addpath(genpath('..'))
model = load('../data/fish_data/fish_X.txt');
scene = load('../data/fish_data/fish_Y.txt');
motion = 'tps';


model(:,3) = 0.2 + model(:,2);
scene(:,3) = 0.1;

%model(1,1)=0.2;
th = pi/4;
tm = [cos(th) -sin(th); sin(th) cos(th)];
%model = model*tm;

[config] = initialize_config(model, scene, motion);
[param, transformed_model, history, config] = gmmreg_L2(config);

figure;

subplot(1,2,1);
DisplayPoints(model,scene,size(model,2)); grid on; %axis off; 
subplot(1,2,2);
DisplayPoints(transformed_model,scene,size(scene,2)); grid on; %axis off
