epsilon = 0.5;
Z = [-1: 0.00001: 1]';

tic;
numel(Z)
X = mexsmooth_poly_C2(Z, epsilon);
DX = mexDsmooth_poly_C2(Z, epsilon);
D2X = mexD2smooth_poly_C2(Z, epsilon);
toc;
%% plot the result
figure;
plot(Z, X, 'b', 'LineWidth', 2);
hold on;
plot(Z, DX, 'r', 'LineWidth', 2);
plot(Z, D2X, 'g', 'LineWidth', 2);
legend('f', 'Df', 'D2f');



