function a = myfunc(inp)

%This is definition using MATLAB's myfunc for solving our initial TDOA
%problem given with the (3-4) in the paper.
a = [0 0];

x = inp(1);
y = inp(2);

a(1) = sqrt((x - 9) * (x - 9) + y * y) - sqrt(x * x + y * y) + (5 - sqrt(52));
a(2) = sqrt((x - 15) * (x - 15) + (y - 7) * (y - 7)) - sqrt(x * x + y * y) + (5 - sqrt(153));