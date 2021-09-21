%this file is for the Monte Carlo approximation for time fractional diffusion Cauchy problem with alpha=1
% du(x,t)/dt = d^2u(x,t)/dx^2, u(x,t=0)=exp(sqrt(2)*x), 
rng(13) 

%define the grid of x (dt), n is the number of stochastic processes (Brownian motion) that we created.
dx=0.01;
x = (-0.5:dx:0.5)';
t = [0.01,0.05,0.1];
n=10000;

% mu is the mean of the normal distribution (the increments of the brownian motion is normal distributed)
length1 = length(x);
length2 = length(t);
mu=0;

% an1, an2, an3 is the analytic solution u(x,t)= exp(sqrt(2)*x+t) at time t=0.01, 0.05, 0.1;
an1 = exp(sqrt(2)*x+t(1));
an2 = exp(sqrt(2)*x+t(2));
an3 = exp(sqrt(2)*x+t(3));

%create 3 cell 
data = cell(length2,1);

%for each cell, sample n data from normal distribution N(0, t)
for i=1:length2
    data{i}= normrnd(mu,sqrt(t(i)),1,n);
end


%create length(x) cell for each specific time 
un1 = cell(length1,1);
for k = 1:length1
    un1{k} = zeros(n,1);
end

%create n brownian motion and apply the monte carlo method
for i=1:length1
    for j=1:n
    un1{i}(j) = q(x(i)+data{1}(j));
    end
end

%take the average
u1 = zeros(length1,1);
for i = 1: length1
    u1(i) = mean(un1{i});
end

for i=1:length1
    enu1=sum(u1(i).*conj(u1(i)).*dx);
end

%generate the figure
figure(1)
plot(x,u1,'d')
hold on 
plot(x,an1)
hold off
xlabel('x');
ylabel('u(x,t=0.01)');
legend("simulation", "analytic");
title('the solution u(x,t) of the diffusion equation at t=0.01');


%the same for other specific time
un2 = cell(length1,1);
for k = 1:length1
    un2{k} = zeros(n,1);
end


for i=1:length1
    for j=1:n
    un2{i}(j) = q(x(i)+data{2}(j));
    end
end

u2 = zeros(length1,1);
for i = 1: length1
    u2(i) = mean(un2{i});
end

for i=1:length1
    enu2=sum(u2(i).*conj(u2(i)).*dx);
end


figure(2)
plot(x,u2,'d')
hold on 
plot(x,an2)
hold off
xlabel('x');
ylabel('u(x,t=0.05)');
legend("simulation", "analytic");
title('the solution u(x,t) of the diffusion equation at t=0.05');


un3 = cell(length1,1);
for k = 1:length1
    un3{k} = zeros(n,1);
end


for i=1:length1
    for j=1:n
    un3{i}(j) = q(x(i)+data{3}(j));
    end
end

u3 = zeros(length1,1);
for i = 1: length1
    u3(i) = mean(un3{i});
end

for i=1:length1
    enu3=sum(u3(i).*conj(u3(i)).*dx);
end


figure(3)
plot(x,u3,'d')
hold on 
plot (x,an3)
hold off
xlabel('x');
ylabel('u(x,t=0.1)');
legend("simulation", "analytic");
title('the solution u(x,t) of the diffusion equation at t=0.1');


figure(4)
plot(x,u1, x,u2,'--',x,u3,'-.')
legend('t=0.001','t=0.05','t=0.1','Location','North')
xlabel('x');
ylabel('u(x,t)');
title('the solution of the fractional diffusion equation u(x,t) at several specific time t');


%compute the error rate by comparing with the analytic values
error_rate1 = zeros(length1,1);
for i=1:length1
    error_rate1(i)=(abs(u1(i)-an1(i)))./an1(i);
end

error_rate2 = zeros(length1,1);
for i=1:length1
    error_rate2(i)=(abs(u2(i)-an2(i)))./an2(i);
end

error_rate3 = zeros(length1,1);
for i=1:length1
    error_rate3(i)=(abs(u3(i)-an3(i)))./an3(i);
end

%compute the maximum error rate
e1 = max(error_rate1);
e2 = max(error_rate2);
e3 = max(error_rate3);


