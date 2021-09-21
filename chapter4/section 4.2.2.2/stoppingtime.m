function stoptime = stoppingtime(alpha,beta,scale,delta,desired_time,n,a)
X = stblrnd(alpha,beta,scale,delta,n,1);
Y = cumsum(X);
B = find(Y>(desired_time-a));
stoptime = Y(B(1));
end