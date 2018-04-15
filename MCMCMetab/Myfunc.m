function [dxdt] = Myfunc(t, Xd,	OptimizeParameters)

k1 = OptimizeParameters(1);
k2 = OptimizeParameters(2);
k3 = OptimizeParameters(3);
k4 = OptimizeParameters(4);
k5 = OptimizeParameters(5);

x = Xd(1);
y = Xd(2);

dxdt(1,1) = k1/(36+k2*y) - k3;
dxdt(2,1) = k4*x - k5;
 
end




