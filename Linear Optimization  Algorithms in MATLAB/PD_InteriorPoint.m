function [xsol, objval] = PD_InteriorPoint(c, A, b)
% MATLAB code for the predictor-corrector method in Ch.6 
% PD_InteriorPoint solves a linear programming in standard form
%                        min  c'*x
%                        s.t. A*x = b
%                             x >= 0
% using the predictor corrector primal-dual path following method
% of Mehrotra.
%
% Inputs:
%  c = n*1 vector, objective coefficients
%  A = m*n matrix with m < n, A is full rank matrix
%  b = m*1 vector, RHS
%
% Outputs:
%  xsol = n*1 vector, final solution
%  objval is  scalar,  final objective value

[m n]=size(A); % number of constraint and variables
e=ones(n,1);
%% Step 0: Initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Obtain an initial interior solution [x(0), pie(0), z(0)]^T such that
% x(0)>0 and z(0)>0. Let k=0 and epsi be some small positive number
% (tolerance). Go to STEP 1.
k=0;                      %counter
epsi=1/10^6;              %tolerance
eta=.95;                  %step length dampening constant
%generate a warm start point
lambda=(A*A')\(2*b); %Lagrange multiplier
x_bar=.5*A'*lambda;  %solve min ||x||+lambda*(b-Ax)
pie_bar=(A*A')\(A*c);%solve min ||A'*pie -c||
z_bar=c-A'*pie_bar;
del_x=max([0; -1.5*min(x_bar)]);
del_z=max([0; -1.5*min(z_bar)]);
del_x_bar=del_x+.5*(x_bar+del_x*e)'*(z_bar+del_z*e)/sum(z_bar+del_z);
del_z_bar=del_z+.5*(x_bar+del_x*e)'*(z_bar+del_z*e)/sum(x_bar+del_x);
x(:,k+1)=x_bar+del_x_bar; %initial x(0), primal variable
pie(:,k+1)=pie_bar;       %initial pie(0), slack variable of dual
z(:,k+1)=z_bar+del_z_bar; %initial z(0), dual variable
obj_pd(:,k+1)=[c'*x(:,k+1); b'*pie(:,k+1)];
Norm(:,k+1)=[norm(A*x(:,k+1)-b); norm(A'*pie(:,k+1)+z(:,k+1)-c); x(:,k+1)'*z(:,k+1);];
while k>=0
    %% Step 1: Affine Scaling Direction Generation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Solve KKT system for affine direction d_affine in the algorithm. GO to STEP 2.
    r_p=A*x(:,k+1)-b;             %primal residuals
    r_d=A'*pie(:,k+1)+z(:,k+1)-c; %dual residuals
    X=diag(x(:,k+1));             %diag(x(k))
    Z=diag(z(:,k+1));             %diag(z(k))
    coeffi_kkt=[zeros(size(A',1), n) A' eye(size(A',1), size(X,2)); ...
        A zeros(m, size(A',2)) zeros(m, size(X,2)); ...
        Z zeros(size(X,1), size(A',2)) X];%coefficient matrix of KKT system
    d_aff=-coeffi_kkt\[r_d; r_p; X*Z*e];  %solve the KKT system
    d_x_aff=d_aff(1:n);                   %affine direction of x(k)
    d_z_aff=d_aff(n+m+1:end);             %affine direction of z(k)
    %% Step 2: Centering Parameter Generation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute alpha_x_affine, alpha_z_affine, y(k), y_affine(k) and let
    % tau(k) = (y_affine(k)/y(k))^3. Solve KKT system for corrector
    % direction d in the algorithm. Go to STEP 3.
    x_temp=x(:,k+1);
    flag_x=find(d_x_aff<0);
    alpha_x_aff = ...
        min([1; min(-x_temp(flag_x)./d_x_aff(flag_x))]);%alpha_x_affine
    z_temp=z(:,k+1);
    flag_z=find(d_z_aff<0);
    alpha_z_aff = ...
        min([1; min(-z_temp(flag_z)./d_z_aff(flag_z))]);%alpha_z_affine
    y(k+1)=x(:,k+1)'*z(:,k+1)/n;   %y(k)
    y_aff(k+1) = ...
        (x(:,k+1)+alpha_x_aff*d_x_aff)'*(z(:,k+1)+alpha_z_aff*d_z_aff)/n;%y_affine(k)
    tau(k+1)=(y_aff(k+1)/y(k+1))^3;%tau(k)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    D_x=diag(d_x_aff);             %D_x(k)
    D_z=diag(d_z_aff);             %D_z(k)
    d = ...
        -coeffi_kkt\[r_d; r_p; X*Z*e+D_x*D_z*e-tau(k+1)*y(k+1)*e];%solve the KKT system
    d_x=d(1:n);                    %d_x
    d_pie=d(n+1:n+m);              %d_pie
    d_z=d(n+m+1:end);              %d_z
    %% Step 3: New Primal and Dual solution Generation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute alpha_x and alpha_z and let x(k+1)=x(k)+alpha_x*d_x,
    % pie(k+1)=pie(k)+alpha_z*d_pie, and z(k+1)=z(k)+alpha_z*d_z.
    % If the stopping criteria is met i.e. ||A*x(k+1) - b||<= epsi,
    % ||A^T*pie(k+1) + z(k+1) - c||<= epsi, and (x(k+1))^T*z(k+1)<= epsi,
    % then STOP. Else k=k+1, go to STEP 1.
    flag_x=find(d_x<0);
    alpha_x_max = ...
        min([1; min(-x_temp(flag_x)./d_x(flag_x))]);% minimum ratio test for x
    alpha_x=min([1; eta*alpha_x_max]); %alpha_x
    flag_z=find(d_z<0);
    alpha_z_max = ...
        min([1; min(-z_temp(flag_z)./d_z(flag_z))]);%minimum ratio test for z
    alpha_z=min([1; eta*alpha_z_max]); %alpha_z
    k=k+1;%update the counter
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    x(:,k+1)=x(:,k)+alpha_x*d_x;               %generate x(k+1)=x(k)+alpha_x*d_x
    pie(:,k+1)=pie(:,k)+alpha_z*d_pie;         %generate pie(k+1)=pie(k)+alpha_z*d_pie
    z(:,k+1)=z(:,k)+alpha_z*d_z;               %generate z(k+1)=z(k)+alpha_z*d_z
    obj_pd(:,k+1)=[c'*x(:,k+1); b'*pie(:,k+1)];%primal and dual objective value
    Norm(:,k+1) = ...
        [norm(A*x(:,k+1)-b); norm(A'*pie(:,k+1)+z(:,k+1)-c); x(:,k+1)'*z(:,k+1);];
    if isempty(find(Norm(:,k+1) >= epsi))%if all norm of residual <= epsi, then optimal STOP.
        disp('problem solved')
        break
    end
end
xsol=x(:,end);         %optimal solution
objval=obj_pd(end,end);%optimal objective value