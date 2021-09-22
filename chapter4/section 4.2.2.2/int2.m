%according to the formula of the Monte Carlo estimator, the variance of the increments of the Brownian motion is changing 
%we define a function int2 to represent t_i^k=(i-1)*k where difference here refers to k.

function intger2 = int2 (x,difference,c)
intger2 = zeros(c,1);
for i=1:x
    intger2(i)= (i-1)*difference;
end
