close all; clear; %clear classes;
clc; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dynamics

SRD = SRDuserinterface;
LinkArray = SRD_get('LinkArray');

%SRD_PrintLinkNames(LinkArray);

mass = sym('m', [3, 1]);

SRD_LinkSet_PropertyByName(LinkArray, 'Shin',  'Mass', mass(1));
SRD_LinkSet_PropertyByName(LinkArray, 'Hip',   'Mass', mass(2));
SRD_LinkSet_PropertyByName(LinkArray, 'Torso', 'Mass', mass(3));


UseCasadi = false;
SymbolicEngine = SRDSymbolicEngine(LinkArray, UseCasadi);
SymbolicEngine.InitializeLinkArray();

SRD_dynamics_derive_JacobiansForLinkArray('SymbolicEngine', SymbolicEngine);

H = SRD_dynamics_derive_JSIM('SymbolicEngine', SymbolicEngine);

[in, dJSIM] = SRD_dynamics_derive_GeneralizedInertialForces_via_dJSIM(...
    'SymbolicEngine', SymbolicEngine, ...
    'JointSpaceInertiaMatrix', H);

g = SRD_dynamics_derive_GeneralizedGravitationalForces(...
    'SymbolicEngine', SymbolicEngine, ...
    'GravitationalConstant', [0; 0; -9.8]);

d = SRD_dynamics_derive_GeneralizedDissipativeForces_uniform(...
    'SymbolicEngine', SymbolicEngine, ...
    'UniformCoefficient', 1);

%NaiveControlMap
T = SRD_dynamics_derive_ControlMap_eye(...
    'SymbolicEngine', SymbolicEngine);






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

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Linearization
% 
% description = SRD_generate_dynamics_linearization(...
%     'SymbolicEngine',                         SymbolicEngine, ...
%     'H', H, ...
%     'c', (in + g + d), ...
%     'T', T, ...
%     'Symbolic_ToOptimizeFunctions',           true, ...
%     'Casadi_cfile_name',                      'g_dynamics_linearization', ...
%     'FunctionName_A',                         'g_linearization_A', ...
%     'FunctionName_B',                         'g_linearization_B', ...
%     'FunctionName_c',                         'g_linearization_c', ...
%     'Path',                                   'Linearization/');
% 
% Handler_dynamics_Linearized_Model = SRD_get_handler__dynamics_linearized_model('description', description);
% SRD_save(Handler_dynamics_Linearized_Model, 'Handler_dynamics_Linearized_Model');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constraints

%%%%%%%%%%%%
%construct constraint
constraint = SymbolicEngine.LinkArray(4).AbsoluteFollower(1);
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
    'dof_robot', SymbolicEngine.dof, ...
    'dof_Constraint', length(constraint));

SRD_save(Handler_Constraints_Model, 'Handler_Constraints_Model');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inverse kinematics

%%%%%%%%%%%%
%construct inverse kinematics task
% rC = SymbolicEngine.GetCoM;
Task = [SymbolicEngine.q(1); SymbolicEngine.q(2); constraint]; 
%%%%%%%%%

description = SRD_generate_second_derivative_Jacobians('SymbolicEngine', SymbolicEngine, ...
    'Task',                                      Task, ...
    'Casadi_cfile_name',                         'g_InverseKinematics', ...
    'Symbolic_ToSimplify',                       true, ...
    'Symbolic_UseParallelizedSimplification',    false, ...
    'Symbolic_ToOptimizeFunctions',              true, ...
    'FunctionName_Task',                         'g_InverseKinematics_Task', ...
    'FunctionName_TaskJacobian',                 'g_InverseKinematics_TaskJacobian', ...
    'FunctionName_TaskJacobian_derivative',      'g_InverseKinematics_TaskJacobian_derivative', ...
    'Path',                                      'InverseKinematics/');

Handler_IK_Model = SRD_get_handler__IK_model('description', description, ...
    'dof_robot', SymbolicEngine.dof, ...
    'dof_Task', length(Task));
SRD_save(Handler_IK_Model, 'Handler_IK_Model');



