%this file is for the Monte Carlo approximation for time fractional diffusion Cauchy problem with alpha to be fractional number alpha=0.5
% d^(0.5)u(x,t)/dt^(0.5) = d^2u(x,t)/dx^2, u(x,t=0)=sin(x), u(x=0,t)=u(x=pi,t)=0, x in [0, pi]

rng(123)

%set the values of alpha, beta, scale ,delta (parameters of stable distributions), n is the number of increments of alpha-stable subordinators 
%(set a large value to ensure the existence of stopping times) 
difference = 0.01;
alpha = 0.5;
beta = 1;
scale = (difference)^(1/alpha);
delta = 0;
n = 100000;
desired_time = [0.01,0.05,0.1];
desired_x = (0:pi/100:pi)';
a = 0;

c = 10000;

gridsize = length(desired_x);

stop_process= cell(c,1);

%generate c alpha-stable subordinators
for i=1:c
    stop_process{i}=cumsum(step(difference,alpha,beta,n));
end

%find T_t^k k=1,...,c
order1 = cell(c,1);
for i=1:c
    order1{i}=find(stop_process{i}>(desired_time(1)-a));
end

stoptime1 = zeros(c,1);
for i=1:c
    stoptime1(i)=stop_process{i}(order1{i}(1));
end
    

%use s function which defined in another file to compute S^k 
bt=cell(gridsize,1);
for i=1:gridsize
    bt{i}=s(desired_x(i),c,difference);
end

btt=cell(gridsize,1);
for i=1:gridsize
    btt{i}=bt{i}*difference;
end

%compute min(S^k, T_t^k) and used it as the variance of the increments of the brownian motion
last=cell(gridsize,1);
for i=1:gridsize
    for j=1:c
        last{i}(j)=min(stoptime1(j),btt{i}(j));
    end
end
    
un1 = cell(gridsize,1);

for k = 1:gridsize
    un1{k} = zeros(c,1);
end

%apply Monte Carlo method to approximate the solutions 
for i=1:gridsize
    for j = 1: c
        un1{i}(j) = g(desired_x(i)+normrnd(0,sqrt(last{i}(j)),1,1));
    end
end

%take the average 
u1 = zeros(gridsize,1);
for i=1:gridsize
    u1(i)= mean(un1{i});
end

%use similar procedures for another desired time.
order2 = cell(c,1);
for i=1:c
    order2{i}=find(stop_process{i}>(desired_time(2)-a));
end

stoptime2 = zeros(c,1);
for i=1:c
    stoptime2(i)=stop_process{i}(order2{i}(1));
end
    

bt1=cell(gridsize,1);
for i=1:gridsize
    bt1{i}=s(desired_x(i),c,difference);
end

btt1=cell(gridsize,1);
for i=1:gridsize
    btt1{i}=bt1{i}*difference;
end

last2=cell(gridsize,1);
for i=1:gridsize
    for j=1:c
        last2{i}(j)=min(stoptime2(j),btt1{i}(j));
    end
end
    
un2 = cell(gridsize,1);

for k = 1:gridsize
    un2{k} = zeros(c,1);
end

for i=1:gridsize
    for j = 1: c
        un2{i}(j) = g(desired_x(i)+normrnd(0,sqrt(last2{i}(j)),1,1));
    end
end

u2 = zeros(gridsize,1);
for i=1:gridsize
    u2(i)= mean(un2{i});
end


%generate the figure of the approximating solution 
figure(1)
plot(desired_x,u1, desired_x,u2,'--')
legend('t=0.01','t=0.05','Location','East')
xlabel('x');
ylabel('u(x,t)');
title('the approximating solution of the fractional diffusion equation u(x,t) at several specific time t');


%compute the analytic values of the solution (based on two-parameter Mittag-leffer function)
k=8000;
vv=zeros(k,1);
for i=1:k
    vv(i)= (((-0.5*(0.01^(0.5)))^(i))/(gamma(0.5*i+1)));
end

ab=sum(vv)+(1/gamma(1));

