%% Performance evaluation

% MATLAB script for finding average time of execution using 100000
% iterations of calling tic and toc of the functions.
numIter = 100000;
options = optimset('Display','off');

xx = zeros(numIter, 1);

for i = 1 : numIter
    tic
    x = lokacijaTDOA(a, b, c, d, e);
    xx(i) = toc;
end

yy = zeros(numIter, 1);

for i = 1 : numIter
    tic
    x = fsolve(@myfunc, [1.5, 2], options);
    yy(i) = toc;
end

disp(['Time - our 2D-TDOA:', num2str(mean(xx))]);
disp(['Time - fsolve in the MATLAB: ', num2str(mean(yy))]);
