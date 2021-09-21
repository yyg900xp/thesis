%this file is for the Monte Carlo approximation for time fractional diffusion Cauchy problem with alpha to be fractional number alpha=0.5
% d^(0.5)u(x,t)/dt^(0.5) = d^2u(x,t)/dx^2+x, u(x,t=0)=sin(x), 


rng(321)
%set the values of alpha, beta (parameters of stable distributions), n is the number of increments of alpha-stable subordinators 
%(set a large value to ensure the existence of stopping times) 
%difference^(1/alpha) is the scale parameter of the stable distributions and c is the number of brownian motion 
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

% use mo2 function which defined in another file to apply the monte carlo approximation u1(x,t=0.01), u2(x,t=0.05).
u1 = zeros(length(x),1);
for i=1:length(x)
    u1(i)=mo2(difference,alpha,beta,n,desired_time(1),x(i),a,c);
end

u2 = zeros(length(x),1);
for i=1:length(x)
    u2(i)=mo2(difference,alpha,beta,n,desired_time(2),x(i),a,c);
end

% generate the figure of the approximating solutions.
figure(1)
plot(x,u1)
xlabel('x');
ylabel('u(x,t=0.01)');
title('the approximating solution of the fractional diffusion equation u(x,t) at t=0.01');




