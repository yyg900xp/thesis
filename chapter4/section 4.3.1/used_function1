function boundtime = s(x,c,difference)
boundary_process= cell(c,1);

for i=1:c
    boundary_process{i}=cumsum(normrnd(0,sqrt(difference),1,5000));
end


bounpx = cell(c,1);
for i=1:c
    bounpx{i}=boundary_process{i}+x;
end

order2 = cell(c,1);
for i=1:c
    order2{i}=find(bounpx{i} <0 | bounpx{i} > pi);
end

order3 = zeros(c,1);
for i=1:c
    order3(i)=order2{i}(1);
end
boundtime=order3;
