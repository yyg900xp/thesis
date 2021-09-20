%This file is for the reproductivity of the paths of the brownian motion B_t

rng(123)
n=100;
dt=0.01;
t = (dt:dt:n*dt)';
t = [0;t];
%define the time interval of the brownian motion t, dt is the time step and n is the number of the increments of the brownian motion.

increments1 = normrnd(0,sqrt(dt),1,n);
increments2 = normrnd(0,sqrt(dt),1,n);
increments3 = normrnd(0,sqrt(dt),1,n);
%create the increments of the brownian motion by sampling data form normal distributions N(0,dt).

path1 = [cumsum(increments1)]';
path1 = [0;path1];
path2 = [cumsum(increments2)]';
path2 = [0;path2];
path3 = [cumsum(increments3)]';
path3 = [0;path3];
%use the cumsum function to create the paths of the brownian motion


figure(1)
plot(t,path1,t,path2,t,path3)
xlabel('time (t)');
ylabel('B_t');
title('Some paths of one dimensional standard Brwonian motion with variance 0.01');
%generate the figures that show the paths of the brownian motion.

