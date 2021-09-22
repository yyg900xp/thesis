%this file is to define a function mo2 to help us obtaining the approximating solution of the time fractional diffusion Cauchy problem 
%by using Monte Carlo method
% the scale and delta refers to the scale and location parameter of the stable distributions 


function montecarlo2 = mo2(difference,alpha,beta,n,desired_time,desired_x,a,c)

scale = (difference)^(1/alpha);
delta = 0;
stop = zeros(c,1);

%compute the stopping time for each alpha-stable subordinator
for i=1:c
    stop(i)= stoppingtime(alpha,beta,scale,delta,desired_time,n,a);
end

%use floor function to compute [T_t/h] (the number of the additions inside the Monte Carlo estimator)
number = floor(stop/difference);

nn = max(number);

%we define a function called int2 to obtain the increasing variance of the increments of the brownian motion.
ele = int2(nn,difference, c);


test = cell(c,1);

for k = 1:c
    test{k} = zeros(number(k),1);
end

%the Monte Carlo approximation for the second part of the Monte Carlo estimator
for i=1:c
    for j = 1: length(test{i})
        test{i}(j) = difference * (normrnd(0,sqrt(ele(j)),1,1)+desired_x);
    end
end

second = zeros(c,1);
for i =1:c
    second(i) = sum(test{i});
end

%g is the initial condition and this is the Monte Carlo approximation for the first part of the Monte Carlo estimator
first = zeros(c,1);
for i = 1:c
    first(i)=g(desired_x+normrnd(0,sqrt(stop(i)),1,1));
end

total = zeros(c,1);
for i = 1:c
    total(i)=first(i)+second(i);
end

%take the average
montecarlo2 = mean(total);

end



        
