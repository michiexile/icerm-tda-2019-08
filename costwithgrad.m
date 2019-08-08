function [f,g] = costwithgrad(x)
% Calculate objective f
%f = 100*(x(2) - x(1)^2)^2 + (1-x(1))^2;
f = sum(abs(fc-delta_Matrix*x) );

if nargout > 1 % gradient required
    g = [-400*(x(2)-x(1)^2)*x(1)-2*(1-x(1));
        200*(x(2)-x(1)^2)];
end
