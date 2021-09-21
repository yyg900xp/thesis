rng(321)
difference = 0.01;
alpha=0.5;
beta=1;
n=1000;
desired_time = [0.01,0.05];
dx=0.1;
x = (0:dx:pi)';
length1=length(x);
a=0;
c=100;

u1 = zeros(length(x),1);
for i=1:length(x)
    u1(i)=mo2(difference,alpha,beta,n,desired_time(1),x(i),a,c);
end

u2 = zeros(length(x),1);
for i=1:length(x)
    u2(i)=mo2(difference,alpha,beta,n,desired_time(2),x(i),a,c);
end

figure(1)
plot(x,u1)
xlabel('x');
ylabel('u(x,t=0.01)');
title('the approximating solution of the fractional diffusion equation u(x,t) at t=0.01');




