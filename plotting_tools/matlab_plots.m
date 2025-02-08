clear;close all;clc

%Import the pendulum data from the FORTRAN Code:
data = readmatrix("ODE_Pendulum_Solution.csv");

%Set the pendulum data from the FORTRAN Code:
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


% %make a graphic animation of the pendulum motion using the FORTRAN data:
% % Define the figure and axis properties
% figure('color','w')
% axis equal
% axis([-1.5 1.5 -1.5 1.5])
% grid on
% hold on
% % Plot the pendulum rod and bob
% rod = line([0, 0], [0, 0], 'Color', 'k', 'LineWidth', 2);
% bob = plot(0, 0, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
% % Loop through the time steps to update the pendulum position
% for i = 1:length(t)
%     % Calculate the pendulum position
%     x = sin(theta(i));
%     y = -cos(theta(i));
%     % Update the rod and bob positions
%     set(rod, 'XData', [0, x], 'YData', [0, y]);
%     set(bob, 'XData', x, 'YData', y);
%     % Pause to create the animation effect
%     pause(0.05);
% end


% % Create a VideoWriter object
% video = VideoWriter('pendulum_animation.mp4','MPEG-4');
% open(video);

% for i = 1:length(t)
%     % Calculate the pendulum position
%     x = sin(theta(i));
%     y = -cos(theta(i));
%     % Update the rod and bob positions
%     set(rod, 'XData', [0, x], 'YData', [0, y]);
%     set(bob, 'XData', x, 'YData', y);
    
%     % Capture the frame and write to the video
%     frame = getframe(gcf);
%     writeVideo(video, frame);
    
%     pause(0.05);
% end

% close(video);

%--------------------------------------------------------------------------
%Import the data from the FORTRAN code for the Damped Harmonic Oscillator:

data_DHO = readmatrix("ODE_Damped_Harmonic_Oscillator_Solution.csv");

%Assign the data to variables:
t_DHO = data_DHO(:,1);
x_DHO = data_DHO(:,2);

%Solve the ODE using OD45 to compare with the imported data:
% Define the initial conditions and time span

span   = [0.0, 25.0];
IC     = [1.0, 0.5];


% Solve the ODE using ode45
[t_ode45_DHO, x_ode45_DHO] = ode45(@Damped_Harmonic_Oscillator, span, IC);
% Plot the results

figure('color','w')
plot(t_DHO, x_DHO, 'LineWidth', 1.5, 'DisplayName','FORTRAN Code');
hold on
plot(t_ode45_DHO, x_ode45_DHO(:,1), 'ro', 'LineWidth', 1.5, 'DisplayName','MATLAB ode45');
xlabel('Time, s')
ylabel('Displacement, m')
legend('Location','best')
grid on




%--------------------------------------------------------------------------------
% ODE Systems: 

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



function xdot = Damped_Harmonic_Oscillator(t, x)
    % Damped_Harmonic_Oscillator
    % ODE system for a damped harmonic oscillator.
    % The equation: x'' + 2*zeta*omega*x' + omega^2*x = 0 is transformed to a
    % first-order system:
    %   x(1) = x      and      x(2) = dx/dt.
    % Thus:
    %   dx/dt    = x(2)
    %   d^2x/dt^2 = -2*zeta*omega*x(2) - omega^2*x(1)
    
    % Define constants for the oscillator
    zeta = 0.3;
    omega = 1.0;
    
    % Ensure that at least two state variables are provided
    if numel(x) < 2
        xdot = [];
        return;
    end
    
    % Preallocate the derivative vector
    xdot = zeros(size(x));
    
    % Compute the derivatives
    xdot(1) = x(2);
    xdot(2) = -2*zeta*omega*x(2) - omega^2*x(1);
end