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
 
C_case1 = [...
     %1 0 0 0 0 0 0 0;
     %0 1 0 0 0 0 0 0;
     0 0 1 0 0 0 0 0;
     0 0 0 1 0 0 0 0;
    % 0 0 0 0 1 0 0 0;
    % 0 0 0 0 0 1 0 0;
     0 0 0 0 0 0 1 0;
     0 0 0 0 0 0 0 1]; 


C = C_case1;
 

ControllerCost = struct('Q', 100*eye(4), 'R', 1);
ObserverCost   = struct('Q', 100*eye(4), 'R', size(C, 1));
System = struct('A', A, 'B', B, 'C', C, 'G', G, 'g', g, 'tol', tol, ...
    'ControllerCost', ControllerCost, 'ObserverCost', ObserverCost, ...
    'x_desired', 0.01*[0;0;randn;randn; 0;0;0;0], ...
    'x_initial', 0.01*[0;0;randn;randn; 0;0;0;0]);
Output = LTI_CLQE(System);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% demo of the node-ability critertion 
N = Output.Matrices.N;
P = (N'*B) * pinv(N'*B);
norm( (eye(size(P)) - P) * (N'*A*N*Output.desired.z_corrected + N'*g) )
norm( (eye(size(P)) - P) * (N'*A*N*Output.desired.z + N'*g) )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % proof that constraint Jacobian does not matter
% size_x = 2;
% size_xi = 2;
% size_l = 1;
% 
% Gs = sym('G', [size_l, size_x]); assume(Gs, 'real');
% Fs = sym('F', [size_x, size_l]); assume(Fs, 'real');
% 
% M =                  [eye(size_x),            zeros(size_x, size_xi),   -Fs;
%                       zeros(size_xi, size_x), eye(size_xi,  size_xi),    zeros(size_xi, size_l);
%                       Gs,                     zeros(size_l, size_xi),    zeros(size_l, size_l)];
% iM = M \ eye(size(M));
% iM11 = iM(1:(size_x+size_xi), 1:(size_x+size_xi))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


G = randn(5, 10);
F = randn(10, 5);
M = [eye(10), -F; G, zeros(5, 5)];
iM = pinv(M);
iM11 = iM(1:10, 1:10)

iM11_v2 = eye(10) - F*pinv(G*F) * G

norm(iM11 - iM11_v2)




