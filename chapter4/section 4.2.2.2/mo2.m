function montecarlo2 = mo2(difference,alpha,beta,n,desired_time,desired_x,a,c)

scale = (difference)^(1/alpha);
delta = 0;
stop = zeros(c,1);

for i=1:c
    stop(i)= stoppingtime(alpha,beta,scale,delta,desired_time,n,a);
end

number = floor(stop/difference);

nn = max(number);

ele = int2(nn,difference, c);


test = cell(c,1);

for k = 1:c
    test{k} = zeros(number(k),1);
end

for i=1:c
    for j = 1: length(test{i})
        test{i}(j) = difference * (normrnd(0,sqrt(ele(j)),1,1)+desired_x);
    end
end

second = zeros(c,1);
for i =1:c
    second(i) = sum(test{i});
end

first = zeros(c,1);
for i = 1:c
    first(i)=g(desired_x+normrnd(0,sqrt(stop(i)),1,1));
end

total = zeros(c,1);
for i = 1:c
    total(i)=first(i)+second(i);
end

montecarlo2 = mean(total);

end



        