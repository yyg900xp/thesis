%the numerical approximation for non-fractional case
rng(13) 

dx=0.01;
x = (-0.5:dx:0.5)';
a=0;
b=0.1;
t = [0.01,0.05,0.1];
n=10000;

length1 = length(x);
length2 = length(t);
mu=0;

an1 = exp(sqrt(2)*x+t(1));
an2 = exp(sqrt(2)*x+t(2));
an3 = exp(sqrt(2)*x+t(3));


data = cell(length2,1);

for i=1:length2
    data{i}= normrnd(mu,sqrt(t(i)),1,n);
end



un1 = cell(length1,1);
for k = 1:length1
    un1{k} = zeros(n,1);
end


for i=1:length1
    for j=1:n
    un1{i}(j) = q(x(i)+data{1}(j));
    end
end

u1 = zeros(length1,1);
for i = 1: length1
    u1(i) = mean(un1{i});
end

for i=1:length1
    enu1=sum(u1(i).*conj(u1(i)).*dx);
end


figure(1)
plot(x,u1,'d')
hold on 
plot(x,an1)
hold off
xlabel('x');
ylabel('u(x,t=0.01)');
legend("simulation", "analytic");
title('the solution u(x,t) of the diffusion equation at t=0.01');

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

e1 = max(error_rate1);
e2 = max(error_rate2);
e3 = max(error_rate3);


