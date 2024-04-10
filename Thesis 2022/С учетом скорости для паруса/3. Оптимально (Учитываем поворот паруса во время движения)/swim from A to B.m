%% clear everything
clc
clear

% Параметры лодки и ветра
boat_k = 10;
boat_k1 = 10;
boat_vb = 5;
boat_m = 1000;
boat_S = 10;

% Начальное положение, скорость и начальное время
boat_pos = 0;
boat_speed = 0;
time = 0;

[boat_x, boat_y, boat_pos ,boat_speed, time] = transfer_between_dots(0,0, -100000,100000,boat_k,boat_k1,boat_vb,boat_m,boat_S, boat_pos, boat_speed,time);
 
%curr_gamma = acot( (0-100)/(0-100) )*180/pi
%curr_gamma = (pi+ acot( (0+100)/(0-100) ) )*180/pi
boat_speed
time
