clear
clc
% MATLAB code for Computational Project Chapter 1 (Section 1.5 pp.39-40)
load ProjectDataCh1
% ProjectDataCh1 is a file that contains monthly per period returns r_it
% in in_return
n = 500; % the number of stocks can be chose
T = 24; % number of time periods
% computing geometric means
for i = 1:n
    mu(i) = (prod(1+in_return(:,i)))^(1/T)-1;
end
R=0:.005: max(mu);          % range of expected return goals R
c = [zeros(n,1); ones(T,1); ones(T,1)]; % compute MAD objective coefficients
Aeq = [];
for t=1:T
    Aeq = cat(1, Aeq, in_return(t,:)-mu);
end
Aeq = [Aeq -eye(T) eye(T);  % constraint coefficients for MAD
    mu zeros(1,2*T);
    ones(1,n) zeros(1,2*T);];
lb = zeros(n+T+T,1);        % lower bound on variables

% computing optimal portfolios over range of return goals R
for a = 1:length(R)
    beq = [zeros(T,1); R(a); 1]; %right hand side coefficients for each R
    [x_MAD(:,a), fval_MAD(a)] = linprog(c, [],[], Aeq,beq, lb,[]);
end

fval_MAD = (1/T)*fval_MAD;  %  minimizing (1/T)w
devi = (pi/2)^.5*fval_MAD;  %  w = sqrt(2/pi)*SD
invest_frac = x_MAD(1:n, :);% optimal portfolio weights for each R

% create figure for optimal portfolios
figure(1)
[xx, yy] = meshgrid(1:n, R);
mesh(xx,yy, invest_frac')
colormap bone
axis([0 500 0 max(mu) 0 1])
xlabel('stocks')
ylabel('expected return R')
zlabel('investment fraction')
title('Portfolio Composition under different R')

% create figure for the efficient frontier of MAD
figure(2)
plot(devi, R, '-k*')
xlabel('volatility \sigma')
ylabel('expected return R')
title('The efficient frontier of MAD')