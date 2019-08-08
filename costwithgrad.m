function [f,g] = costwithgrad(x)
% Calculate objective f
delta_Matrix=dlmread('D:/icerm-tda-2019-08/coboundary.txt');
%
fc=dlmread('D:/icerm-tda-2019-08/l2cocycle.txt');
z0=zeros(size(delta_Matrix,1),1);
delta_Matrix=transpose(delta_Matrix);

%f = 100*(x(2) - x(1)^2)^2 + (1-x(1))^2;
f = sum(abs(fc-delta_Matrix*x) );

if nargout > 1 % gradient required
    g = [-400*(x(2)-x(1)^2)*x(1)-2*(1-x(1));
        200*(x(2)-x(1)^2)];
end
