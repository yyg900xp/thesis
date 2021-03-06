%the code here is to create the histograms of the stable distributions under different chosen alpha.
%the algorithm that used here is CMS method (Chamber's method) which described in the thesis.
%the chosen parameters are alpha=(0.3, 0.5, 0.8), beta = 1, sigma = 1^alpha, mu = 0.
%the number of sampled data is 1000.

rng(123)
n=1000;
alpha1=0.3;
alpha2=0.5;
alpha3=0.8;
beta=1;
dt=1;
sigma1=dt^(alpha1);
sigma2=dt^(alpha2);
sigma3=dt^(alpha3);
mu=0;
% set the values of the parameters of the stable distributions.


%apply the CMS method
%give the values of zeta and epsilon which defined in the theis.
zeta1 = - beta * tan(pi*alpha1/2);
epsilon1 = (1/alpha1)*atan(-zeta1);
zeta2 = - beta * tan(pi*alpha2/2);
epsilon2 = (1/alpha2)*atan(-zeta2);
zeta3 = - beta * tan(pi*alpha3/2);
epsilon3 = (1/alpha3)*atan(-zeta3);


U = unifrnd(-0.5*pi,0.5*pi,1,n);
W = exprnd(1,n,1);
%sample n independent data from uniform distribution (-0.5 pi, 0.5 pi) and exponential distribution (mean 1).


%create n data that are stable distorbuted S(alpha1,beta,1,0)
X1=zeros(n,1);
for i=1:n
    X1(i)=((1+zeta1^2)^(1/(2*alpha1)))*(sin(alpha1*(U(i)+epsilon1))/((cos(U(i)))^(1/alpha1)))...
        *((cos(U(i)-alpha1*(U(i)+epsilon1))/(W(i)))^((1-alpha1)/alpha1));
end 

%create n data that are stable distorbuted S(alpha,beta,sigma1,mu)
Y1=zeros(n,1);
for i=1:n
    Y1(i)=sigma1*X1(i)+mu; 
end



Y4 = Y1(find(Y1 < 20));


%create n data that are stable distorbuted S(alpha2,beta,1,0)
X2=zeros(n,1);
for i=1:n
    X2(i)=((1+zeta2^2)^(1/(2*alpha2)))*(sin(alpha2*(U(i)+epsilon2))/((cos(U(i)))^(1/alpha2)))...
        *((cos(U(i)-alpha2*(U(i)+epsilon2))/(W(i)))^((1-alpha2)/alpha2));
end 

%create n data that are stable distorbuted S(alpha,beta,sigma2,mu)
Y2=zeros(n,1);
for i=1:n
    Y2(i)=sigma2*X2(i)+mu; 
end


Y5 = Y2(find(Y2 < 20));

%create n data that are stable distorbuted S(alpha3,beta,1,0)
X3=zeros(n,1);
for i=1:n
    X3(i)=((1+zeta3^2)^(1/(2*alpha3)))*(sin(alpha3*(U(i)+epsilon1))/((cos(U(i)))^(1/alpha3)))...
        *((cos(U(i)-alpha3*(U(i)+epsilon3))/(W(i)))^((1-alpha3)/alpha3));
end 

%create n data that are stable distorbuted S(alpha,beta,sigma3,mu)
Y3=zeros(n,1);
for i=1:n
    Y3(i)=sigma3*X3(i)+mu; 
end

Y6 = Y3(find(Y3 < 20));

% create n data for each chosen alpha case that are stable distributed by CMS method.
% the formula are follow the definition which are given in the thesis.
% we use the find function to just look at the behaviors of the small sampled data of the histogrmas since the
% stable distributions will generate some large data which affect the histograms (looks like a pulse function).

figure(1)
histogram(Y4,50);
xlabel('Y');
ylabel('Number of data');
ylim([0,300]);
title('Histogram of the sampled data by CMS method');


figure(2)
histogram(Y5,50);
xlabel('Y');
ylabel('Number of data');
ylim([0,300]);
title('Histogram of the sampled data by CMS method');

figure(3)
histogram(Y6,70);
xlabel('Y');
ylabel('Number of data');
ylim([0,300]);
title('Histogram of the sampled data by CMS method');

% create the histograms of the stable distributions under different chosen alpha.

