%this file is for showing the paths of the alpha-stable subordinator under different chosen alpha.


%define the values of the parameters alpha=(a1,a2,a3), beta, sigma=dt^alpha, mu=delta=0
%the increments of the alpha-stable subordinator is sampled from S(alpha,beta,sigma,delta)
%the number of the increments is n
%the time interval is t and the difference of time is dt.
a1 = 0.1;
a2 = 0.5;
a3 = 0.8;
beta=1;
dt=0.01;
n=1000;
delta=0;
t = (dt:dt:n*dt)';

%use step function which defined in another file to sample the data that are stable distributed.
increments1=step(dt,a1,beta,n);
increments2=step(dt,a2,beta,n);
increments3=step(dt,a3,beta,n);

%use cumsum function to create the path of the alpha-stable subordinator.
path1=cumsum(increments1);
path2=cumsum(increments2);
path3=cumsum(increments3);

%generate the corresponding figures of the path of the alpha-stable subordinator.
figure(1)
plot(t,path1,'d','MarkerSize',3)
xlabel('time t');
ylabel('X_t');
title('the path of the alpha-stable Levy motion X_t with alpha=0.1');


figure(2)
plot(t,path2,'d','MarkerSize',3)
xlabel('time t');
ylabel('X_t');
title('the path of the alpha-stable Levy motion X_t with alpha=0.5');


figure(3)
plot(t,path3,'d','MarkerSize',3)
xlabel('time t');
ylabel('X_t');
title('the path of the alpha-stable Levy motion X_t with alpha=0.8');
