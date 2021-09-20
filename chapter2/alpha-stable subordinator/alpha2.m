
a1 = 0.1;
a2 = 0.5;
a3 = 0.8;
beta=1;
dt=0.01;
n=1000;
delta=0;
t = (dt:dt:n*dt)';

increments1=step(dt,a1,beta,n);
increments2=step(dt,a2,beta,n);
increments3=step(dt,a3,beta,n);

path1=cumsum(increments1);
path2=cumsum(increments2);
path3=cumsum(increments3);


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
