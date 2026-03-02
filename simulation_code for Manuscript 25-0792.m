clear; clc; close all;

% System matrices
S = [0 1; -1 0];
M = [0 1; -2 -3];
N = [0; 1];

% Event-triggered parameters
xi = 0.01;
delta = 1e-3;
beta = 5;
alpha = 1;

% Simulation time
t_start = 0;
t_end   = 100;      % adjustable as needed

% State indices (for readability)
ind_v1 = 1;
ind_v2 = 2;
ind_z1 = 3;
ind_z2 = 4;
ind_y  = 5;
ind_eta1 = 6;
ind_eta2 = 7;
ind_psi1 = 8;
ind_psi2 = 9;
ind_K   = 10;
nstates = 10;

%% Initial conditions (according to homework.pdf point 3)
% v(0) = [2; -0.8]
v1_0 = 2;
v2_0 = -0.8;
% z1(0), z2(0), y(0)
z1_0 = -0.6;
z2_0 = -0.8;
y_0  = -0.2;
% eta(0)
eta1_0 = -0.7;
eta2_0 = -0.5;
% psi(0)
psi1_0 = 0.6;
psi2_0 = -0.1;
% K(0)
K_0 = 2.5;

% Combined initial state vector
x0 = [v1_0; v2_0; z1_0; z2_0; y_0; eta1_0; eta2_0; psi1_0; psi2_0; K_0];

% Parameters at the initial event instant (t=0 as the first event)
e_prev    = y_0;
rho_prev  = e_prev^4 + 1;
K_prev    = K_0;
psi_prev  = [psi1_0, psi2_0];   % row vector
eta_prev  = [eta1_0; eta2_0];    % column vector

% Initial control input (held constant during the first event interval)
u_const = -K_prev * rho_prev * e_prev + psi_prev * eta_prev;

% Record all event instants (including t=0)
event_times = 0;

%% Main loop: piecewise integration until simulation ends
% Store results from all time segments
all_t = [];
all_x = [];
current_t = t_start;
current_x = x0;

% Set ODE options (event function to be defined later)
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8, 'Events', ...
    @(t,x) event_condition(t, x, xi, delta, beta, alpha, ...
    K_prev, e_prev, rho_prev, psi_prev, eta_prev, ...
    ind_y, ind_K, ind_eta1, ind_eta2, ind_psi1, ind_psi2));

while current_t < t_end
    % Define dynamics for the current interval (capturing u_const and indices)
    f = @(t,x) system_dynamics(t, x, u_const, M, N, ...
        ind_v1, ind_v2, ind_z1, ind_z2, ind_y, ...
        ind_eta1, ind_eta2, ind_psi1, ind_psi2, ind_K);
    
    % Integrate from current_t to t_end
    [T, X, te, ye] = ode45(f, [current_t, t_end], current_x, options);
    
    % Append current results to the total data (avoid duplicating the first point)
    if isempty(all_t)
        all_t = T;
        all_x = X;
    else
        all_t = [all_t; T(2:end)];
        all_x = [all_x; X(2:end, :)];
    end
    
    % Check if an event is triggered
    if ~isempty(te)
        % Record event instant
        event_times = [event_times; te(end)];
        
        % Update previous values (using the state ye at the event end)
        current_x = ye(end, :)';   % new starting state
        e_prev    = current_x(ind_y);
        rho_prev  = e_prev^4 + 1;
        K_prev    = current_x(ind_K);
        psi_prev  = [current_x(ind_psi1), current_x(ind_psi2)];
        eta_prev  = [current_x(ind_eta1); current_x(ind_eta2)];
        
        % Update constant control input
        u_const = -K_prev * rho_prev * e_prev + psi_prev * eta_prev;
        
        % Update current time
        current_t = te(end);
        
        % Reset the event function (capturing new previous values)
        options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8, 'Events', ...
            @(t,x) event_condition(t, x, xi, delta, beta, alpha, ...
            K_prev, e_prev, rho_prev, psi_prev, eta_prev, ...
            ind_y, ind_K, ind_eta1, ind_eta2, ind_psi1, ind_psi2));
    else
        % No event triggered, simulation ends
        break;
    end
end

% Extract time vector and states
t = all_t;
v1 = all_x(:, ind_v1);
v2 = all_x(:, ind_v2);
z1 = all_x(:, ind_z1);
z2 = all_x(:, ind_z2);
y  = all_x(:, ind_y);
eta1 = all_x(:, ind_eta1);
eta2 = all_x(:, ind_eta2);
psi1 = all_x(:, ind_psi1);
psi2 = all_x(:, ind_psi2);
K    = all_x(:, ind_K);

% Tracking error e = y
e = y;

%% Plotting
% Figure 1: Tracking error e(t)
figure(1);
plot(t, e, 'b', 'LineWidth', 1.0);
xlabel('Time {\itt} (sec)');ylabel('$\varphi_{1}$');
set(gcf, 'Position', [100, 100, 400, 230])

% Figure 2: States z1(t) and z2(t)
figure(2);
plot(t, z1+1, 'b', 'LineWidth', 1.0); hold on;
plot(t, z2+1, 'r', 'LineWidth', 1.0);
h=legend('$\varphi_{2}$','$\varphi_{3}$');
set(h,'Interpreter','Latex');
set(gcf, 'Position', [100, 100, 400, 230])
xlabel('Time {\itt} (sec)');

% Figure 3: Event-triggered mechanism 
figure(3);
% Red envelope: beta * exp(-alpha * t)
t_smooth = linspace(0, 5, 1000);
beta_curve = beta * exp(-alpha * t_smooth);
plot(t_smooth, beta_curve, 'r', 'LineWidth', 1.5); hold on;
% Blue sawtooth: xi*(||vc||^2 - delta*e^2)
plot(fig3_t, fig3_xi, 'b', 'LineWidth', 1.0);
h = legend('$\beta \mathrm{e}^{-\alpha t}$', '$\xi(\|\tilde{v}_c\|^2 - \delta e^2)$');
set(h, 'Interpreter', 'Latex');
xlabel('Time {\itt} (sec)'); ylabel('Event-triggered intervals');
xlim([0 5]); ylim([0 inf]);
set(gcf, 'Position', [100, 100, 400, 230]);
hold off;

% Figure 4: K(t)
figure(4);
plot(t, K, 'b', 'LineWidth', 1.0);
xlabel('Time {\itt} (sec)');ylabel('\it{K(t)}');
set(gcf, 'Position', [100, 100, 400, 230])

% Figure 5: Two components of \hat{\Psi}(t)
figure(5);
plot(t, psi1, 'b', 'LineWidth', 1.0); hold on;
plot(t, psi2, 'r', 'LineWidth', 1.0);
h=legend('$\hat{\Psi}_{1}$','$\hat{\Psi}_{2}$');
set(h,'Interpreter','Latex');
set(gcf, 'Position', [100, 100, 400, 230])
xlabel('Time {\itt} (sec)');ylabel('$\hat{\Psi}(t)$');


%% Auxiliary function definitions

% System dynamics function (includes leader, follower, controller)
function dxdt = system_dynamics(t, x, u_const, M, N, ...
    ind_v1, ind_v2, ind_z1, ind_z2, ind_y, ...
    ind_eta1, ind_eta2, ind_psi1, ind_psi2, ind_K)
    dxdt = zeros(length(x), 1);
    
    % Current state
    e = x(ind_y);
    eta = [x(ind_eta1); x(ind_eta2)];
    
    % Leader dynamics (1)
    dxdt(ind_v1) = x(ind_v2);
    dxdt(ind_v2) = -x(ind_v1);
    
    % Follower dynamics (2)
    dxdt(ind_z1) = -0.5 * x(ind_z1) + 0.5 * e;
    dxdt(ind_z2) = -0.5 * x(ind_z2) - 0.5 * e;
    dxdt(ind_y)  = e - (1/3)*e^3 - x(ind_z1) + x(ind_z2) + x(ind_v1) + u_const;
    
    % \eta dynamics (3)
    deta = M * eta + N * u_const;
    dxdt(ind_eta1) = deta(1);
    dxdt(ind_eta2) = deta(2);
    
    % \hat{\Psi} dynamics (3)
    dxdt(ind_psi1) = -eta(1) * e;
    dxdt(ind_psi2) = -eta(2) * e;
    
    % K dynamics (3)
    rho = e^4 + 1;
    dxdt(ind_K) = rho * e^2;
end

% Event-triggered condition function (5) and (6)
function [value, isterminal, direction] = event_condition(t, x, xi, delta, beta, alpha, ...
    K_prev, e_prev, rho_prev, psi_prev, eta_prev, ...
    ind_y, ind_K, ind_eta1, ind_eta2, ind_psi1, ind_psi2)
    % Current state
    e = x(ind_y);
    K = x(ind_K);
    psi = [x(ind_psi1), x(ind_psi2)];
    eta = [x(ind_eta1); x(ind_eta2)];
    rho = e^4 + 1;
    
    % Compute the two components of \tilde{v}_c (6) 
    vc1 = psi_prev * eta_prev - psi * eta;          % scalar
    vc2 = -K_prev * rho_prev * e_prev + K * rho * e; % scalar
    norm_vc_sq = vc1^2 + vc2^2;
    
    % Event condition value: triggered when this value crosses from negative to positive
    value = xi * (norm_vc_sq - delta * e^2) - beta * exp(-alpha * t);
    
    isterminal = 1;   % Stop integration when event is detected
    direction = 1;    % Detect only rising direction
end
