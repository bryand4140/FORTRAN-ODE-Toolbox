clear;close all;clc

data = readmatrix("ODE_Pendulum_Solution.csv");



t = data(:,1);
theta = data(:,2);  


%Solve the ODE using OD45 to compare with the imported data:
% Define the initial conditions and time span
IC = [pi/4, 0.0];
t_span = [0.0, 5.0];

% Solve the ODE using ode45
[t_ode45, x_ode45] = ode45(@Pendulum_ODE, t_span, IC);

% Plot the results
figure('color','w')
plot(t, theta, 'LineWidth', 1.5, 'DisplayName','FORTRAN Code');
hold on
plot(t_ode45, x_ode45(:,1), 'ro', 'LineWidth', 1.5, 'DisplayName','MATLAB ode45');
xlabel('Time, s')
ylabel('Angle, rad')
legend('Location','best')
grid on


%make a graphic video of the pendulum motion using the FORTRAN data:
% Define the figure and axis properties
figure('color','w')
axis equal
axis([-1.5 1.5 -1.5 1.5])
grid on
hold on
% Plot the pendulum rod and bob
rod = line([0, 0], [0, 0], 'Color', 'k', 'LineWidth', 2);
bob = plot(0, 0, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
% Loop through the time steps to update the pendulum position
for i = 1:length(t)
    % Calculate the pendulum position
    x = sin(theta(i));
    y = -cos(theta(i));
    % Update the rod and bob positions
    set(rod, 'XData', [0, x], 'YData', [0, y]);
    set(bob, 'XData', x, 'YData', y);
    % Pause to create the animation effect
    pause(0.05);
end


% Create a VideoWriter object
video = VideoWriter('pendulum_animation.mp4','MPEG-4');
open(video);

for i = 1:length(t)
    % Calculate the pendulum position
    x = sin(theta(i));
    y = -cos(theta(i));
    % Update the rod and bob positions
    set(rod, 'XData', [0, x], 'YData', [0, y]);
    set(bob, 'XData', x, 'YData', y);
    
    % Capture the frame and write to the video
    frame = getframe(gcf);
    writeVideo(video, frame);
    
    pause(0.05);
end

close(video);



function xdot = Pendulum_ODE(t, x)
    % ODE system for a simple pendulum (small or large amplitude).
    % The pendulum equation: theta'' + (g/L)*sin(theta) = 0 is transformed to 
    % a first-order system:
    %   x(1) = theta      and      x(2) = dtheta/dt.
    % Thus, we have:
    %   dtheta/dt = x(2)
    %   d^2theta/dt^2 = - (g/L) * sin(x(1))
    
    % Define constants for gravity and pendulum length.
    g = 9.81;
    L = 1.0;
    
    % Ensure that at least two state variables are available.
    if numel(x) < 2
        error('Input vector x must have at least 2 elements.');
    end
    
    % Preallocate the derivative vector.
    xdot = zeros(2, 1);
    
    xdot(1) = x(2);
    xdot(2) = - (g / L) * sin(x(1));
end