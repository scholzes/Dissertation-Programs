function [y] = sinc(x)

y = sin(x)./x;
y(isnan(y)) = 1;

end

