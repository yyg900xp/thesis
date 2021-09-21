rng(123)
difference = 0.0001;
alpha = 0.5;
beta = 1;
n = 100000;
desired_time = [0.01,0.05,0.1];
desired_x = (0:pi/1000:pi)';
a = 0;

c = 10000;

gridsize = length(desired_x);

process= cell(c,1);

for i=1:c
    process{i}=cumsum(step(difference,alpha,beta,n));
end

order1 = cell(c,1);
for i=1:c
    order1{i}=find(process{i}>(desired_time(1)-a));
end

stoptime1 = zeros(c,1);
for i=1:c
    stoptime1(i)=process{i}(order1{i}(1));
end
    
un1 = cell(gridsize,1);

for k = 1:gridsize
    un1{k} = zeros(c,1);
end

for i=1:gridsize
    for j = 1: c
        un1{i}(j) = g(desired_x(i)+normrnd(0,sqrt(stoptime1(j)),1,1));
    end
end

u1 = zeros(gridsize,1);
for i=1:gridsize
    u1(i)= mean(un1{i});
end


order2 = cell(c,1);
for i=1:c
    order2{i}=find(process{i}>(desired_time(2)-a));
end

stoptime2 = zeros(c,1);
for i=1:c
    stoptime2(i)=process{i}(order2{i}(1));
end

un2 = cell(gridsize,1);

for k = 1:gridsize
    un2{k} = zeros(c,1);
end

for i=1:gridsize
    for j = 1: c
        un2{i}(j) = g(desired_x(i)+normrnd(0,sqrt(stoptime2(j)),1,1));
    end
end

u2 = zeros(gridsize,1);
for i=1:gridsize
    u2(i)= mean(un2{i});
end


order3 = cell(c,1);
for i=1:c
    order3{i}=find(process{i}>(desired_time(3)-a));
end

stoptime3 = zeros(c,1);
for i=1:c
    stoptime3(i)=process{i}(order3{i}(1));
end
    
un3 = cell(gridsize,1);

for k = 1:gridsize
    un3{k} = zeros(c,1);
end

for i=1:gridsize
    for j = 1: c
        un3{i}(j) = g(desired_x(i)+normrnd(0,sqrt(stoptime3(j)),1,1));
    end
end

u3 = zeros(gridsize,1);
for i=1:gridsize
    u3(i)= mean(un3{i});
end

figure(1)
plot(desired_x,u1, desired_x,u2,'--',desired_x,u3,'-.')
legend('t=0.01','t=0.05','t=0.1','Location','East')
xlabel('x');
ylabel('u_N^k(x,t)');
title('the approximating solutions u_N^k(x,t) of the fractional diffusion equation at several specific time')




