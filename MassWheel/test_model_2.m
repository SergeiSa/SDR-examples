close all; clear classes; clc;

tol = 10^(-5);

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

H = Handler_dynamics_generalized_coordinates_model.get_joint_space_inertia_matrix(q);
iH = pinv(H);
% c = Handler_dynamics_generalized_coordinates_model.get_bias_vector(q0, v0);

A = Handler_dynamics_Linearized_Model.get_A(q, v, u, iH);
B = Handler_dynamics_Linearized_Model.get_B(q, v,    iH);
% C = [eye(n), zeros(n, n)];
%C = [eye(n+1), zeros(n+1, n-1)];
% % C = [0 0 1 0 0 0 0 0;
% %      0 0 0 1 0 0 0 0;
% %      0 0 0 0 0 0 1 0;
% %      0 0 0 0 0 0 0 1];
C = [...
     %1 0 0 0 0 0 0 0;
     %0 1 0 0 0 0 0 0;
     0 0 1 0 0 0 0 0;
     0 0 0 1 0 0 0 0;
    % 0 0 0 0 1 0 0 0;
    % 0 0 0 0 0 1 0 0;
     0 0 0 0 0 0 1 0;
     0 0 0 0 0 0 0 1];
% C = randn(3, 2*n);


F  = Handler_Constraints_Model.get_Jacobian(q);
dF = Handler_Constraints_Model.get_Jacobian_derivative(q);
G = [zeros(k, n), F; F, dF];
N = null(G); R = orth(G'); %E = [N, R];

mm = size(C, 1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%D1*dx/dt + D2*x = 0
D1 = [eye(n), zeros(n)];
D2 = [zeros(n), eye(n)];

GG = [G, zeros(2*k, 2*n); D1, D2];
GG = svd_suit(GG, tol);

V = GG.null; V = V((2*n+1):end, :);

Ra = orth(V*pinv(V)*R);
%Rn = orth((eye(2*n) - V*pinv(V))*R);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NA = svd_suit(N'*A);
Re = NA.row_space;

temp = svd_suit( (Re*Re' * V*pinv(V) * R), tol);
Rae = temp.orth;

nn_1 = size(N,  2);
nn_2 = size(Rae, 2);
nn = nn_1+nn_2;

E = [N, Rae];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Kz    = lqr(N'*A*N, N'*B, eye(nn_1), eye(m));
Kzeta = pinv(N'*B) * N'*A*Rae;
K = [Kz, Kzeta];

A1 = [N'*A*N, N'*A*Rae; zeros(nn_2, nn_1), zeros(nn_2, nn_2)];
A2 = [N'*B*K; zeros(nn_2, nn)];

L = lqr(A1', E'*C', eye(nn), eye(mm));
L = L';

Com1 = [N'*A*N-N'*B*Kz,   N'*B*K;
        zeros(nn, nn_1), (A1 - L*C*E)];
eig(Com1)    

Com2 = [N'*A*N, -N'*B*K;
        L*C*N,  (A1 - A2 - L*C*E)];
eig(Com2)    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% cvx_solver Gurobi_2
% cvx_begin sdp
%     
%     variable L2(2*k, mm)
%     
%     minimize(0)
% 
%     subject to
%         L2 * C * N == zeros(k, (2*n - 2*k));
%         L2 * C * R == semidefinite
% cvx_end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% x = 0.1*randn(2*n, 1); z = N'*x; zeta = R'*x; xi_est = 0.1*randn(2*n, 1); Y0 = [z; xi_est]; 
x = [0;0;randn;randn; 0;0;0;0]; z = N'*x; zeta = Rae'*x; xi_est = 0.1*randn(nn, 1); Y0 = [z; xi_est]; 

v = [N'*A*Rae; L*C*Rae] * zeta;
ODEFUN = get_fnc(Com2, v);

[TOUT,YOUT] = ode45(ODEFUN,[0 100], Y0);
figure('Color', 'w')
subplot(2, 2, 1)
plot(TOUT,YOUT(:, 1:nn_1))
subplot(2, 2, 2)
plot(TOUT,YOUT(:, (1+nn_1):(2*nn_1) ))
subplot(2, 2, 3)

zeta = YOUT(:, (1+2*nn_1):(2*nn_1+nn_2) );
if ~isempty(zeta)
    plot(TOUT,YOUT(:, (1+2*nn_1):(2*nn_1+nn_2) ))
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fnc = get_fnc(Com2, v)

    fnc = @my_fnc;
    function dy = my_fnc(~, y)
        dy = Com2 * y + v;
    end

end

