rng(123)
n=100;
dt=0.01;
t = (dt:dt:n*dt)';
t = [0;t];
increments1 = normrnd(0,sqrt(dt),1,n);
increments2 = normrnd(0,sqrt(dt),1,n);
increments3 = normrnd(0,sqrt(dt),1,n);

path1 = [cumsum(increments1)]';
path1 = [0;path1];
path2 = [cumsum(increments2)]';
path2 = [0;path2];
path3 = [cumsum(increments3)]';
path3 = [0;path3];


figure(1)
plot(t,path1,t,path2,t,path3)
xlabel('time (t)');
ylabel('B_t');
title('Some paths of one dimensional standard Brwonian motion with variance 0.01');


