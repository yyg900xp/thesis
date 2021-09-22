function intger2 = int2 (x,difference,c)
intger2 = zeros(c,1);
for i=1:x
    intger2(i)= (i-1)*difference;
end
