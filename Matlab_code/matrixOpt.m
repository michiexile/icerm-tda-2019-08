delta_Matrix=dlmread('D:/icerm-tda-2019-08/coboundary.txt');
%
f=dlmread('D:/icerm-tda-2019-08/l2cocycle.txt');
z0=zeros(size(delta_Matrix,1),1);

delta_Matrix=transpose(delta_Matrix);

%costfunc is L^2.
costfunc2 = @(x)transpose(f-delta_Matrix*x)*(f-delta_Matrix*x);
%costfunc is L^1.
costfunc1 = @(x)sum(abs(f-delta_Matrix*x) );
%costfunc is L^p.
p=4;
costfuncp = @(x)sum( (f-delta_Matrix*x).^p );
%Set up lambda.
lambda=0.4
%costfunc is (1-lambda)*L^1+lambda*L^2
costfunc12 = @(x)(1-lambda)*costfunc1(x)+lambda*costfunc2(x);
%costfunc is L^1+lambda*L^p
costfunc1p = @(x)costfunc1(x)+lambda*costfuncp(x);

%options = optimoptions('fmincon','Display','iter');
%z_min = fminsearch(costfunc12,z0);%derivative free method, very slow.
%z_min = fmincon(costfunc12,z0,transpose(z0),10);sum(z_m); % constrained
%minimizer, fast bu tinaccurate.

%options = optimoptions('fminunc','Algorithm','trust-region','SpecifyObjectiveGradient',true);
options = optimoptions(@fminunc,'Display','iter','Algorithm','quasi-newton','MaxIterations',20);

%x0 = [-1,2];
%fun = @costwithgrad;
%z_min = fminunc(costfunc1,z0,options)
z_min = fminunc(costfunc12,z0,options)


%Show the minimum and the argmin.
z_min
costf_min = costfunc12(z_min)

%Display the function values.
scatter(point_cloud(:,3), point_cloud(:,4), '.')
axis equal
