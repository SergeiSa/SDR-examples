close all; clear; clc;

tol = 10^(-5);

% [0, 0, phi_1, phi_2]
q = [0; 0; 0; 0];
v = zeros(size(q));
x = [q; v];
u = 0;

Handler_dynamics_generalized_coordinates_model = SRD_get('Handler_dynamics_generalized_coordinates_model');
Handler_dynamics_Linearized_Model              = SRD_get('Handler_dynamics_Linearized_Model_c');
Handler_Constraints_Model                      = SRD_get('Handler_Constraints_Model');


k = Handler_Constraints_Model.dof_Constraint;
dof = Handler_dynamics_generalized_coordinates_model.dof_configuration_space_robot;
m = Handler_dynamics_generalized_coordinates_model.dof_control;

H = Handler_dynamics_generalized_coordinates_model.get_joint_space_inertia_matrix(q);
c = Handler_dynamics_generalized_coordinates_model.get_bias_vector(q, v);
T = Handler_dynamics_generalized_coordinates_model.get_control_map();

F  = Handler_Constraints_Model.get_Jacobian(q);
dF = Handler_Constraints_Model.get_Jacobian_derivative(q);

iH = pinv(H);
f0 = [v; iH*(T*u - c)];

M = [H, -F';
     F, zeros(k,k)];
iM = pinv(M);

A = Handler_dynamics_Linearized_Model.get_A(q, v, u, iM);
B = Handler_dynamics_Linearized_Model.get_B(q, v,    iM);
% A = Handler_dynamics_Linearized_Model.get_A(q, v, u, iH);
% B = Handler_dynamics_Linearized_Model.get_B(q, v,    iH);
g = f0 - A * x - B * u;

G = [zeros(k, dof), F; F, dF];



C_case0 = [...
    % 1 0 0 0 0 0 0 0;
    % 0 1 0 0 0 0 0 0;
    % 0 0 1 0 0 0 0 0;
     0 0 0 1 0 0 0 0;
    % 0 0 0 0 1 0 0 0;
    % 0 0 0 0 0 1 0 0;
    % 0 0 0 0 0 0 1 0;
    % 0 0 0 0 0 0 0 1...
    ]; 
 
C_case1 = [...
     %1 0 0 0 0 0 0 0;
     %0 1 0 0 0 0 0 0;
     0 0 1 0 0 0 0 0;
     0 0 0 1 0 0 0 0;
    % 0 0 0 0 1 0 0 0;
    % 0 0 0 0 0 1 0 0;
     0 0 0 0 0 0 1 0;
     0 0 0 0 0 0 0 1]; 
 
C_case2 = [...
     1 0 0 0 0 0 0 0;
     0 1 0 0 0 0 0 0;
     0 0 1 0 0 0 0 0;
     0 0 0 1 0 0 0 0;
    % 0 0 0 0 1 0 0 0;
    % 0 0 0 0 0 1 0 0;
     0 0 0 0 0 0 1 0;
     0 0 0 0 0 0 0 1];
 
C_case3 = [...
     1 0 0 0 0 0 0 0;
     0 1 0 0 0 0 0 0;
     0 0 1 0 0 0 0 0;
     0 0 0 1 0 0 0 0;
    % 0 0 0 0 1 0 0 0;
     0 0 0 0 0 1 0 0;
    % 0 0 0 0 0 0 1 0;
     0 0 0 0 0 0 0 1]; 

C_case4 = eye(2*dof); 

C = C_case1;
 

ControllerCost = struct('Q', 100*eye(4), 'R', 1);
ObserverCost   = struct('Q', 100*eye(4), 'R', size(C, 1));
System = struct('A', A, 'B', B, 'C', C, 'G', G, 'g', g, 'tol', tol, ...
    'ControllerCost', ControllerCost, 'ObserverCost', ObserverCost, ...
    'x_desired', 0.01*[0;0;randn;randn; 0;0;0;0], ...
    'x_initial', 0.01*[0;0;randn;randn; 0;0;0;0]);
Output = LTI_CLQE(System);


Time = 15;
[TOUT,YOUT] = ode45(Output.closed_loop.z_xi.ode_fnc, [0 Time], Output.closed_loop.z_xi.Y0);

z        = YOUT(:, 1:Output.sizes.size_z);
z_est    = YOUT(:, (1+Output.sizes.size_z):(2*Output.sizes.size_z));
zeta_est = YOUT(:, (1+2*Output.sizes.size_z):(2*Output.sizes.size_z+Output.sizes.size_zeta) );
x_est = (Output.Matrices.N *z_est' + Output.Matrices.R_used*zeta_est')';

figure('Color', 'w')
subplot(2, 2, 1)
plot(TOUT, z, 'LineWidth', 1.5); hold on; title('z');
plot([TOUT(1); TOUT(end)], [Output.desired.z_corrected'; Output.desired.z_corrected'], '--', 'LineWidth', 0.8);

subplot(2, 2, 2)
plot(TOUT, z_est); title('z est')
subplot(2, 2, 3)
if ~isempty(zeta_est)
    plot(TOUT,zeta_est); title('zeta est')
end

subplot(2, 2, 4)

plot(TOUT, x_est, 'LineWidth', 1.5); hold on; title('x est');
plot([TOUT(1); TOUT(end)], [Output.desired.x_corrected'; Output.desired.x_corrected'], '--', 'LineWidth', 0.8);


drawnow;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[TOUT,YOUT] = ode45(Output.closed_loop.x_xi.ode_fnc, [0 Time], Output.closed_loop.x_xi.Y0);

x        = YOUT(:, 1:Output.sizes.size_x);
z_est    = YOUT(:, (1+Output.sizes.size_x):(Output.sizes.size_x+Output.sizes.size_z));
zeta_est = YOUT(:, (1+Output.sizes.size_x+Output.sizes.size_z):(Output.sizes.size_x+Output.sizes.size_z+Output.sizes.size_zeta) );
z_calc = (Output.Matrices.N'*x')';

figure('Color', 'w')
subplot(3, 2, 1)
plot(TOUT, x, 'LineWidth', 1.5); hold on; title('x');
plot([TOUT(1); TOUT(end)], [Output.desired.x_corrected'; Output.desired.x_corrected'], '--', 'LineWidth', 0.8);

subplot(3, 2, 2)
plot(TOUT, z_est); title('z est')
subplot(3, 2, 3)

if ~isempty(zeta_est)
    plot(TOUT, zeta_est); title('zeta est')
end

subplot(3, 2, 4)
plot(TOUT, z_calc); hold on; title('z calc')
plot([TOUT(1); TOUT(end)], [Output.desired.z_corrected'; Output.desired.z_corrected'], '--', 'LineWidth', 0.8);
subplot(3, 2, 5)
plot(TOUT, z-z_calc); title('z-z calc')

drawnow;



