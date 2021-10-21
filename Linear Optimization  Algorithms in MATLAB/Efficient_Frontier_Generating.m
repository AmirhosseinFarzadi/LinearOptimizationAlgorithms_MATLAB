clear
clc
%%%%% Three asset MVO problem in Example 7.3 in Chapter 7 %%%%%
n=3;
%%%%% Data for MVO problem %%%%%
mu=[9.73 6.57 5.37]/100; % expected returns of assets
Q=[.02553 .00327 .00019; %covariance matrix
    .00327 .01340 -.00027;
    .00019 -.00027 .00125];
goal_R=[5.5:.5:9.5]/100; % expected return goals range from 5.5% to 9.5%
for a=1:length(goal_R)
    c=zeros(n,1);
    A=-mu;
    b=-goal_R(a);
    Aeq=[ones(1,n);];
    beq=[1;];
    %%%%% quadratic optimization call %%%%%
    [x(a,:), fval(a,1)] = quadprog(Q, c, A,b, Aeq,beq, [],[]);
    std_devi(a,1)=(2*fval(a,1))^.5; %standard deviation = (x'*Q*x)^.5
end
%%%%% efficient frontier plot %%%%%
plot(std_devi, goal_R, '-k*')
xlabel('volatility \sigma')
ylabel('expected return goal R')
title('The efficient frontier of MVO')