function [y] = MyObjectiveFunction(x)
k = numel(x);
Sum = 0;
for i = 1:k
	xi = x(i);
	Sum = Sum + (xi)^2;
end
y = Sum;
end
