%this file is for computing the stopping time T_t =inf{s>0, tau_s -a >t } where tau_s is the alpha-stable subordinator 
%the stblrnd is a function which can generate data that are stable distributed S(alpha,beta,scale,delta).
%Y is the path of the alpha-stable subordinator.

function stoptime = stoppingtime(alpha,beta,scale,delta,desired_time,n,a)
X = stblrnd(alpha,beta,scale,delta,n,1);
Y = cumsum(X);
B = find(Y>(desired_time-a));
stoptime = Y(B(1));
end
