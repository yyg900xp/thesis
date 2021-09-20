%This file is used to generate a data that is stable distributed S(alpha, beta,difference^(1/alpha),mu).
%The algorithm is based on CMS method which shown in the thesis.
%the n refers to the number of the data that we want to sample.

function Y = step(difference,alpha,beta,n)
sigma=difference^(1/alpha);
mu=0;
zeta = - beta * tan(pi*alpha/2);
epsilon = (1/alpha)*atan(-zeta);


U = unifrnd(-0.5*pi,0.5*pi,1,n);
W = exprnd(1,n,1);

X=zeros(n,1);
for i=1:n
    X(i)=((1+zeta^2)^(1/(2*alpha)))*(sin(alpha*(U(i)+epsilon))/((cos(U(i)))^(1/alpha)))...
        *((cos(U(i)-alpha*(U(i)+epsilon))/(W(i)))^((1-alpha)/alpha));
end 

Y=zeros(n,1);
for i=1:n
    Y(i)=sigma*X(i)+mu; 
end
end

