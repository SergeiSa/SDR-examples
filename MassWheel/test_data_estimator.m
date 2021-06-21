close all; clear; clc;

tol = 10^(-5);

% [0, 0, phi_1, phi_2]
q = [0; 0; pi; 0];
v = zeros(size(q));
x = [q; v];
u = 0;

Handler_dynamics_generalized_coordinates_model = SRD_get('Handler_dynamics_generalized_coordinates_model');
Handler_dynamics_Linearized_Model              = SRD_get('Handler_dynamics_Linearized_Model_c');
Handler_Constraints_Model                      = SRD_get('Handler_Constraints_Model');


k = Handler_Constraints_Model.dof_Constraint;
n = Handler_dynamics_generalized_coordinates_model.dof_configuration_space_robot;
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


C_case0 = [...
    % 1 0 0 0 0 0 0 0;
    % 0 1 0 0 0 0 0 0;
    % 0 0 1 0 0 0 0 0;
     0 0 0 1 0 0 0 0;
    % 0 0 0 0 1 0 0 0;
    % 0 0 0 0 0 1 0 0;
    % 0 0 0 0 0 0 1 0;
     0 0 0 0 0 0 0 1 ...
    ]; 

 
G = [zeros(k, n), F; F, dF];
N = null(G); R = orth(G'); %E = [N, R];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%D1*dx/dt + D2*x = 0
D1 = [eye(n), zeros(n)];
D2 = [zeros(n), eye(n)];

GG = [G, zeros(2*k, 2*n); D1, D2];
GG = svd_suit(GG, tol);

V = GG.null; V = V((2*n+1):end, :);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

temp2 = svd_suit( ( pinv(N'*A)*(N'*A) * V*pinv(V) * R), tol);
Rae = temp2.orth;

our_R = Rae; C = C_case0;

E = [N, our_R];

nn_1 = size(N,  2);
nn_2 = size(our_R, 2);
nn = nn_1+nn_2;
mm = size(C, 1);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Kz    = lqr(N'*A*N, N'*B, 100*eye(nn_1), 1*eye(m));
Kzeta = pinv(N'*B) * N'*A*our_R;
K = [Kz, Kzeta];

N1 = [N, zeros(2*n, nn_2)];

A1 = N1'*A*E;
B1 = N1'*B;

% A1 = [N'*A*N,            N'*A*our_R; 
%       zeros(nn_2, nn_1), zeros(nn_2, nn_2)];
% B1 = [N'*B; 
%       zeros(nn_2, m)];

L = lqr(A1', E'*C', diag([100, 10, 100, 10]), 100*diag([10, 10]));
L = L';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data = readtable('data_static.csv');
data = table2array(data);

data_time = data(:, 1);
data_state = data(:, 2:5);
% data(:, 6) is empty
dataU     = data(:, 8);

data_IMU = data(:, 10:11);
data_obs = data(:, 12:15);
%[data.Var11, data.Var12, data.Var13, data.Var14];



data_phi1 = data_IMU(:, 1);
data_phi2 = data_state(:, 2);
data_w1   = data_IMU(:, 2);
data_w2   = data_state(:, 4);
dataY = [data_phi2, data_w2];

data_state_mesurements = [data_phi1, data_phi2, data_w1, data_w2];

max_time = max(data_time);

current_time = get_current_time(data_time);
current_y = get_current_y(data_time, dataY);
current_u = get_current_u(data_time, dataU);



Est_dynamics = @(t, xi) L*current_y(t) + (A1 - L*C*E)*xi - B1*current_u(t) + N1'*g;
Y0 = [N'*x, our_R'*x];

[TOUT,YOUT] = ode45(Est_dynamics,[0 max_time], Y0);

z_est       = YOUT(:, 1:nn_1);
zeta_est    = YOUT(:, (nn_1+1):(nn_1+nn_2));
x_est = (N*z_est' + our_R*zeta_est')';

figure('Color', 'w')
subplot(2, 2, 1)
plot(TOUT, z_est, 'LineWidth', 1.5); hold on; title('z est');
subplot(2, 2, 2)
if ~isempty(zeta_est)
plot(TOUT, zeta_est, 'LineWidth', 1.5); hold on; title('zeta est');
end
subplot(2, 2, 3)
plot(TOUT, x_est, 'LineWidth', 1.5); hold on; title('x est');
drawnow

figure('Color', 'w')
subplot(2, 2, 1)
plot(TOUT, x_est(:, 3:4), 'LineWidth', 1.5); hold on; title('q1, q2 est');
subplot(2, 2, 2)
plot(data_time, data_state_mesurements(:, 1:2), 'LineWidth', 1.5); hold on; title('q1, q2 inac mes');
subplot(2, 2, 3)
plot(data_time, data_obs(:, [1, 2]), 'LineWidth', 1.5); hold on; title('q1, q2 inac obs');
drawnow


figure('Color', 'w')
subplot(2, 2, 1)
plot(TOUT, x_est(:, 7:8), 'LineWidth', 1.5); hold on; title('v1, v2 est');
subplot(2, 2, 2)
plot(data_time, data_state_mesurements(:, 3:4), 'LineWidth', 1.5); hold on; title('v1, v2 inac mes');
subplot(2, 2, 3)
plot(data_time, data_obs(:, [3, 4]), 'LineWidth', 1.5); hold on; title('v1, v2 inac obs');
drawnow


figure('Color', 'w')
subplot(2, 2, 1)
plot(TOUT, x_est(:, 3), 'LineWidth', 1.5); 
hold on; title('q1');
plot(data_time, data_obs(:, 1), '--', 'LineWidth', 1.5);
plot(data_time, data_phi1, ':', 'LineWidth', 2.5);
legend('est', 'obs', 'imu')
drawnow
















function current_time = get_current_time(data_time)
current_time = @ct;
    function time = ct(t)
        [~, index] = min(abs(data_time - t));
        time = data_time(index);
    end
end
function current_y = get_current_y(data_time, dataY)
current_y = @cy;
    function Y = cy(t)
        [~, index] = min(abs(data_time - t));
        Y = reshape(dataY(index, :), [], 1);
    end
end
function current_u = get_current_u(data_time, dataU)
current_u = @cu;
    function U = cu(t)
        [~, index] = min(abs(data_time - t));
        U = dataU(index, :);
    end
end
