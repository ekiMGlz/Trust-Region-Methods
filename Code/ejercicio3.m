%% Function definition, initial point and optimal point

beale = @(x)  (1.5 - x(1, :) + x(1, :) .* x(2, :)).^2    ...
            + (2.25 - x(1, :) + x(1, :) .* x(2, :).^2).^2 ...
            + (2.625 - x(1, :) + x(1, :) .* x(2, :).^3).^2;

x_0 = [-3; -3];
x_opt = [3; 0.5];

%% Calculate the trajectories

% Total number of iterations for each method to converge on x_opt. The
% number of iterations is greater than the ones previously found, since a
% constant trust region radius was used.
Cauchy_iters = 57;
Dogleg_iters = 17;

% Initialize empty (num_iters)x2 matrices to store both trajectories
Cauchy_trajectory = zeros(Cauchy_iters, 2);
Dogleg_trajectory = zeros(Dogleg_iters, 2);

% Cauchy Trajectory
x_k = x_0;
for i = 1:Cauchy_iters
    Cauchy_trajectory(i,:) = x_k;
    [x_k, ~, ~] = mRC1(beale, x_k, 1);
end

% Dogleg Trajectory
x_k = x_0;
for i = 1:Dogleg_iters
    Dogleg_trajectory(i,:) = x_k;
    [x_k, ~,~] = mRC2(beale, x_k, 1);
end

%% Plotting

hold on
grid on
axis([-4.5, 4.5, -4.5, 4.5])
daspect([1, 1, 1])

% Plot the contour of the beale function 
high_levels = [beale(x_0), 7000, 6000, 5000, 4000, 3000, 2000, 1000, 100, 50];
low_levels = 0:1:10;

colormap parula
fcontour(@(x,y) beale([x;y]), 'LevelList', high_levels);
colorbar

% Plot the trajectories
l1 = plot(Cauchy_trajectory(:, 1), Cauchy_trajectory(:, 2), '-x', 'Color', '#0072BD');
l2 = plot(Dogleg_trajectory(:, 1), Dogleg_trajectory(:, 2), '-x', 'Color', '#D95319');

% Plot initial and final points
plot(x_0(1), x_0(2),'o', 'Color', '	#A2142F', 'MarkerSize', 10);
plot(x_opt(1), x_opt(2),'kh', 'MarkerSize', 10);
legend([l1, l2], {'Trayectoria usando Punto de Cauchy', 'Trayectoria usando Dogleg'});

hold off
uiwait

% Plot a zoomed in version of the trajectory with mRC1
hold on
grid on
axis([0, 3.5, -2, 1.5])
daspect([1, 1, 1])

fcontour(@(x,y) beale([x;y]), 'LevelList', low_levels);
colorbar

% Plot trajectories
plot(Cauchy_trajectory(:, 1), Cauchy_trajectory(:, 2), '-x', 'Color', '#0072BD');

% Plot optimal point
plot(x_opt(1), x_opt(2),'kh', 'MarkerSize', 10);

hold off
uiwait

% Plot a zoomed in version of the trajectory with mRC2
hold on
grid on
axis([0, 3.5, -2, 1.5])
daspect([1, 1, 1])

fcontour(@(x,y) beale([x;y]), 'LevelList', low_levels);
colorbar

% Plot trajectories
plot(Dogleg_trajectory(:, 1), Dogleg_trajectory(:, 2), '-x', 'Color', '#D95319');

% Plot optimal point
plot(x_opt(1), x_opt(2),'kh', 'MarkerSize', 10);

hold off