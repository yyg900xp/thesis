%This file is to use Crank Nicolson method to approximate the solutions of the time fractional diffusion Cauchy problem.
%d^(0.5)u(x,t)/dt^(0.5) = d^2u(x,t)/dx^2, u(x,t=0)=sin(x), x in [0, pi], u(x=0,t)=u(x=pi,t)=0

n=101;
t=linspace(0,1,n);% set time and distance scales
x=(0:pi/100:pi);
m=length(x);
dx=x(2)-x(1);dt=t(2)-t(1); %Increments
s=dt/(dx^2);%Useful for the solution.
s1 = 0.5+ (0.25)*dt+0.5*s;
s2 = -0.25*s;
s3 = -0.25*s;

u=zeros(n,m); %set up solution matrix
A=zeros(m,m);
A(1,1) = 1;

for j=2:m-1
  A(j,j-1) = s2;
  A(j,j) = s1;
  A(j,j+1) = s3;
end
A(m,m) = 1;

%Add in intial condition:
u(1,:)=sin(x);
v=zeros(m,1);
%Solve the system
for i=2:n-1
  %Construct the RHS for the solution
  for j=2:m-1
      v(j)=0.25*s*u(i-1,j+1)+(0.5-(0.25)*dt-0.5*s)*u(i-1,j)+0.25*s*u(i-1,j-1)+0.5*sin(x(j))*dt;
  end
  %Solve for the new time step
  w=A\v;
  u(i,:)=w;
end

%generate the figure of the approximating solutions.
figure(1)
plot(x,u(2,:),x,u(6,:),'--')
legend('t=0.01','t=0.05','Location','East')
xlabel('x');
ylabel('u(x,t)');
title('the approximating solution of the fractional diffusion equation u(x,t) at several specific time t');


% compute the corresponding analytic solutions
k=8000;
vv=zeros(k,1);
for i=1:k
    vv(i)= (((-0.5*(0.01^(0.5)))^(i))/(gamma(0.5*i+1)));
end

ab=sum(vv)+(1/gamma(1));
