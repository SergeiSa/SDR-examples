
close all; clear; %clear classes;
clc; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dynamics

LinkArray = SRD_get('LinkArray');

SymbolicEngine = SRDSymbolicEngine('LinkArray', LinkArray, 'Casadi', false);
SymbolicEngine.InitializeLinkArray();

SRD_dynamics_derive_JacobiansForLinkArray('SymbolicEngine', SymbolicEngine);

H = SRD_dynamics_derive_JSIM('SymbolicEngine', SymbolicEngine);

[in, dH] = SRD_dynamics_derive_GeneralizedInertialForces_via_dH(...
    'SymbolicEngine', SymbolicEngine, ...
    'JointSpaceInertiaMatrix', H);

g = SRD_dynamics_derive_GeneralizedGravitationalForces(...
    'SymbolicEngine', SymbolicEngine, ...
    'GravitationalConstant', [0; 0; -9.8]);

% d = SRD_dynamics_derive_GeneralizedDissipativeForces_uniform(...
%     'SymbolicEngine', SymbolicEngine, ...
%     'UniformCoefficient', 1);
% 
% T = sym([0;0;0;1]);

d = sym([0;0;0;0.00015]);
% d = sym([1;1;0;0.00015]);

T = sym([0;0;0;0.061]);

description = SRD_generate_dynamics_generalized_coordinates_model(...
    'SymbolicEngine', SymbolicEngine, ...
    'H', H, ...
    'c', (in + g + d), ...
    'T', T, ...
    'Symbolic_ToOptimizeFunctions', true, ...
    'Casadi_cfile_name', 'g_dynamics_generalized_coordinates', ...
    'FunctionName_H', 'g_dynamics_H', ...
    'FunctionName_c', 'g_dynamics_c', ...
    'FunctionName_T', 'g_dynamics_T', ...
    'Path', 'Dynamics/');

Handler_dynamics_generalized_coordinates_model = SRD_get_handler__dynamics_generalized_coordinates_model('description', description);
SRD_save(Handler_dynamics_generalized_coordinates_model, 'Handler_dynamics_generalized_coordinates_model');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Linearization

description = SRD_generate_dynamics_linearization(...
    'SymbolicEngine',                         SymbolicEngine, ...
    'H', H, ...
    'c', (in + g + d), ...
    'T', T, ...
    'Symbolic_ToOptimizeFunctions',           true, ...
    'Casadi_cfile_name',                      'g_dynamics_linearization', ...
    'FunctionName_A',                         'g_linearization_A', ...
    'FunctionName_B',                         'g_linearization_B', ...
    'FunctionName_c',                         'g_linearization_c', ...
    'Path',                                   'Linearization/');

Handler_dynamics_Linearized_Model = SRD_get_handler__dynamics_linearized_model('description', description);
SRD_save(Handler_dynamics_Linearized_Model, 'Handler_dynamics_Linearized_Model');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constraints

%%%%%%%%%%%%
%construct constraint
constraint = [SymbolicEngine.q(1); SymbolicEngine.q(2)];
%%%%%%%%%

description = SRD_generate_second_derivative_Jacobians('SymbolicEngine', SymbolicEngine, ...
    'Task',                                   constraint, ...
    'Casadi_cfile_name',                     'g_Constraints', ...
    'Symbolic_ToSimplify',                    true, ...
    'Symbolic_UseParallelizedSimplification', false, ...
    'Symbolic_ToOptimizeFunctions',           true, ...
    'FunctionName_Task',                     'g_Constraint', ...
    'FunctionName_TaskJacobian',             'g_Constraint_Jacobian', ...
    'FunctionName_TaskJacobian_derivative',  'g_Constraint_Jacobian_derivative', ...
    'Path',                                  'Constraints/');


Handler_Constraints_Model = SRD_get_handler__Constraints_model('description', description, ...
    'dof_robot', SymbolicEngine.dof);
SRD_save(Handler_Constraints_Model, 'Handler_Constraints_Model');


F = jacobian(constraint, SymbolicEngine.q);
dF = jacobian(F(:), SymbolicEngine.q) * SymbolicEngine.v;
dF = reshape(dF, size(F));

description = SRD_generate_dynamics_linearization_c(...
    'SymbolicEngine',                         SymbolicEngine, ...
    'H', H, ...
    'c', (in + g + d), ...
    'T', T, ...
    'F', F, ...
    'dF', dF, ...
    'Symbolic_ToOptimizeFunctions',           true, ...
    'Casadi_cfile_name',                      'g_dynamics_linearization', ...
    'FunctionName_A',                         'g_linearization_A', ...
    'FunctionName_B',                         'g_linearization_B', ...
    'FunctionName_c',                         'g_linearization_c', ...
    'Path',                                   'Linearization_c/');

Handler_dynamics_Linearized_Model_c = SRD_get_handler__dynamics_linearized_model('description', description);
SRD_save(Handler_dynamics_Linearized_Model_c, 'Handler_dynamics_Linearized_Model_c');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



tol = 10^(-5);

% [0, 0, phi_1, phi_2]
q = [0; 0; pi/6; 0];
v = zeros(size(q));
x = [q; v];
u = 0;

Handler_dynamics_generalized_coordinates_model = SRD_get('Handler_dynamics_generalized_coordinates_model');
Handler_dynamics_Linearized_Model              = SRD_get('Handler_dynamics_Linearized_Model');
Handler_dynamics_Linearized_Model_c             = SRD_get('Handler_dynamics_Linearized_Model_c');
Handler_Constraints_Model                      = SRD_get('Handler_Constraints_Model');


k = Handler_Constraints_Model.dof_Constraint;
n = Handler_dynamics_generalized_coordinates_model.dof_configuration_space_robot;
m = Handler_dynamics_generalized_coordinates_model.dof_control;

H = Handler_dynamics_generalized_coordinates_model.get_joint_space_inertia_matrix(q);
c = Handler_dynamics_generalized_coordinates_model.get_bias_vector(q, v);
T = [0;0;0;1];
iH = pinv(H);
f0 = [v; iH*(T*u - c)];

A = Handler_dynamics_Linearized_Model.get_A(q, v, u, iH);
B = Handler_dynamics_Linearized_Model.get_B(q, v,    iH);

F  = Handler_Constraints_Model.get_Jacobian(q);
dF = Handler_Constraints_Model.get_Jacobian_derivative(q);
G = [zeros(k, n), F; F, dF];
N = null(G); R = orth(G'); %E = [N, R];

M = [H, -F';
     F, zeros(k,k)];
iM = pinv(M)

Ac = Handler_dynamics_Linearized_Model_c.get_A(q, v, u, iM);
Bc = Handler_dynamics_Linearized_Model_c.get_B(q, v,    iM);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ac
Bc

N'*Ac*N
Kz  = lqr(N'*Ac*N, N'*Bc, eye(4), eye(1))
eig(N'*Ac*N - N'*Bc*Kz)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A
N'*A*N

get_a = getgetter(Handler_dynamics_generalized_coordinates_model);

delta = 0.001;
deltaq = eye(4)*delta;

a  = get_a(q, v, u);
a1 = get_a(q+deltaq(:, 1), v, u);
a2 = get_a(q+deltaq(:, 2), v, u);
a3 = get_a(q+deltaq(:, 3), v, u);
a4 = get_a(q+deltaq(:, 4), v, u);

a5 = get_a(q, v+deltaq(:, 1), u);
a6 = get_a(q, v+deltaq(:, 2), u);
a7 = get_a(q, v+deltaq(:, 3), u);
a8 = get_a(q, v+deltaq(:, 4), u);

AA = ...
[a - a1, a - a2, a - a3, a - a4, ...
 a - a5, a - a6, a - a7, a - a8]/delta;

AA = [zeros(4), eye(4); AA]
N'*AA*N


get_ac = getgetter2(Handler_dynamics_generalized_coordinates_model, Handler_Constraints_Model);

a  = get_ac(q, v, u);
a1 = get_ac(q+deltaq(:, 1), v, u);
a2 = get_ac(q+deltaq(:, 2), v, u);
a3 = get_ac(q+deltaq(:, 3), v, u);
a4 = get_ac(q+deltaq(:, 4), v, u);

a5 = get_ac(q, v+deltaq(:, 1), u);
a6 = get_ac(q, v+deltaq(:, 2), u);
a7 = get_ac(q, v+deltaq(:, 3), u);
a8 = get_ac(q, v+deltaq(:, 4), u);

AA = ...
[a - a1, a - a2, a - a3, a - a4, ...
 a - a5, a - a6, a - a7, a - a8]/delta;

AA = [zeros(4), eye(4); AA]
N'*AA*N

b = get_ac(q, v, u+delta);

BB = ...
(a - b)/delta;
BB = [zeros(4, 1); BB]


Kz  = lqr(N'*AA*N, N'*BB, eye(4), eye(1))
eig(N'*AA*N - N'*B*Kz)

function getter = getgetter(Handler_dynamics_generalized_coordinates_model)
getter = @get_a;
    function a = get_a(q, v, u)
        H = Handler_dynamics_generalized_coordinates_model.get_joint_space_inertia_matrix(q);
        c = Handler_dynamics_generalized_coordinates_model.get_bias_vector(q, v);
        T = [0;0;0;1];
        a = pinv(H)*(T*u - c);
    end
end

function getter = getgetter2(Handler_dynamics_generalized_coordinates_model, Handler_Constraints_Model)
getter = @get_a;
    function a = get_a(q, v, u)
        H = Handler_dynamics_generalized_coordinates_model.get_joint_space_inertia_matrix(q);
        c = Handler_dynamics_generalized_coordinates_model.get_bias_vector(q, v);
        T = [0;0;0;1];
        F = Handler_Constraints_Model.get_Jacobian(q);
        dF = Handler_Constraints_Model.get_Jacobian_derivative(q);
        k = Handler_Constraints_Model.dof_Constraint;
        n = Handler_dynamics_generalized_coordinates_model.dof_configuration_space_robot;
        
        M = [H, -F';
             F, zeros(k, k)];
         
        a = [eye(n), zeros(n, k)] * pinv(M) * [(T*u - c); -dF*v];
    end
end
