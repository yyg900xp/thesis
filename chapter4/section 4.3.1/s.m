%this file is to define a function f for creating c Brownian motion and computing min(S,T_t) as the variance
%of the increments of the Brownian motion where S=inf{s > 0,B_s^x not in [m,n] } [m,n] is the boundary of the solution u(x,t)

function boundtime = s(x,c,difference)
boundary_process= cell(c,1);

%create c standard brownian motion 
for i=1:c
    boundary_process{i}=cumsum(normrnd(0,sqrt(difference),1,5000));
end


%create c path of the brownian motion with different starting point x 
bounpx = cell(c,1);
for i=1:c
    bounpx{i}=boundary_process{i}+x;
end

%compute the first time that the brownian motion escape from the domain [m,n]
order2 = cell(c,1);
for i=1:c
    order2{i}=find(bounpx{i} <0 | bounpx{i} > pi);
end

order3 = zeros(c,1);
for i=1:c
    order3(i)=order2{i}(1);
end
boundtime=order3;
