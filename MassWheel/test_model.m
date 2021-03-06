close all; clear classes; clc;

% [0, 0, phi_1, phi_2]
q = reshape(SRD_get('InitialPosition'), [], 1);
v = zeros(size(q));
u = 0;

Handler_dynamics_generalized_coordinates_model = SRD_get('Handler_dynamics_generalized_coordinates_model');
Handler_dynamics_Linearized_Model              = SRD_get('Handler_dynamics_Linearized_Model');
Handler_Constraints_Model                      = SRD_get('Handler_Constraints_Model');


k = Handler_Constraints_Model.dof_Constraint;
n = Handler_dynamics_generalized_coordinates_model.dof_configuration_space_robot;
m = Handler_dynamics_generalized_coordinates_model.dof_control;
nn_1 = 2*n - 2*k;
nn_2 = 2*k;

H = Handler_dynamics_generalized_coordinates_model.get_joint_space_inertia_matrix(q);
iH = pinv(H);
% c = Handler_dynamics_generalized_coordinates_model.get_bias_vector(q0, v0);

A = Handler_dynamics_Linearized_Model.get_A(q, v, u, iH);
B = Handler_dynamics_Linearized_Model.get_B(q, v,    iH);
% C = [eye(n), zeros(n, n)];
%C = [eye(n+1), zeros(n+1, n-1)];
C = [0 0 1 0 0 0 0 0;
     0 0 0 1 0 0 0 0;
     0 0 0 0 0 0 1 0;
     0 0 0 0 0 0 0 1];
C = [...
   %  1 0 0 0 0 0 0 0;
     0 1 0 0 0 0 0 0;
     0 0 1 0 0 0 0 0;
     0 0 0 1 0 0 0 0;
     0 0 0 0 1 0 0 0;
     0 0 0 0 0 1 0 0;
     0 0 0 0 0 0 1 0;
     0 0 0 0 0 0 0 1];
% C = randn(3, 2*n);


F  = Handler_Constraints_Model.get_Jacobian(q);
dF = Handler_Constraints_Model.get_Jacobian_derivative(q);
G = [zeros(k, n), F; F, dF];
N = null(G); R = orth(G'); E = [N, R];
I = eye(2*n); N = I(:, [3,4,7,8]); R = I(:, [1,2,5,6]); E = [N, R];

mm = size(C, 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Kz    = lqr(N'*A*N, N'*B, eye(2*n - 2*k), eye(m));
Kzeta = pinv(N'*B) * N'*A*R;
K = [Kz, Kzeta];

A1 = [N'*A*N, N'*A*R; zeros(nn_2, nn_1), zeros(nn_2, nn_2)];
A2 = [N'*B*K; zeros(nn_2, 2*n)];

L     = lqr(A1', E'*C', eye(2*n), eye(mm));
L = L';

Com1 = [N'*A*N-N'*B*Kz,   N'*B*K;
        zeros(2*n, nn_1), (A1 - L*C*E)];
eig(Com1)    

Com2 = [N'*A*N, -N'*B*K;
        L*C*N,  (A1 - A2 - L*C*E)];
eig(Com2)    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% x = 0.1*randn(2*n, 1); z = N'*x; zeta = R'*x; xi_est = 0.1*randn(2*n, 1); Y0 = [z; xi_est]; 
x = [0;0;randn;randn; 0;0;0;0]; z = N'*x; zeta = R'*x; xi_est = 0.1*randn(2*n, 1); Y0 = [z; xi_est]; 

v = [N'*A*R; L*C*R] * zeta;
ODEFUN = get_fnc(Com2, v);

[TOUT,YOUT] = ode45(ODEFUN,[0 100], Y0);
figure('Color', 'w')
subplot(2, 2, 1)
plot(TOUT,YOUT(:, 1:nn_1))
subplot(2, 2, 2)
plot(TOUT,YOUT(:, (1+nn_1):(1+2*nn_1) ))
subplot(2, 2, 3)
plot(TOUT,YOUT(:, (1+2*nn_1):end ))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fnc = get_fnc(Com2, v)

    fnc = @my_fnc;
    function dy = my_fnc(~, y)
        dy = Com2 * y + v;
    end

end

