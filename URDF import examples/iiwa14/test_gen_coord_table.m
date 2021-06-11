close all; clear classes; clc;

InitialPosition = SRD_get('InitialPosition');

Handler_State = SRD_get_handler__state('InitialPosition', InitialPosition, ...
    'InitialVelocity', zeros(size(InitialPosition)));

Handler_IK_Solution = SRD_get('Handler_IK_Solution');

Handler_Time = SRDHandler_Time();

Handler_Desired_State = SRD_get_handler__desired_state(...
    'Handler_ControlInput', Handler_IK_Solution, ...
    'Handler_Time',   Handler_Time);

% %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Handler_dynamics_generalized_coordinates_model = SRD_get('Handler_dynamics_generalized_coordinates_model');
Handler_dynamics_Linearized_Model = SRD_get('Handler_dynamics_Linearized_Model');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tf = Handler_IK_Solution.TimeExpiration;
tf = 0.1;
time_table = 0:0.005:tf;

[H_table, iH_table, T_table, c_table, q_table, v_table, a_table, u_table] = ...
    SRD_GenCoordinatesModel_GenerateTable('Handler_dynamics_generalized_coordinates_model', Handler_dynamics_generalized_coordinates_model, ...
    'Handler_IK_Solution', Handler_IK_Solution, ...
    'TimeTable', time_table);


% %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = Handler_dynamics_generalized_coordinates_model.dof_configuration_space_robot;

Q = 1000*eye(2 * n);
R = 0.01*eye(Handler_dynamics_generalized_coordinates_model.dof_control);
Count = size(A_table, 3);
K_table = SRD_LQR_GenerateTable(A_table, B_table, repmat(Q, [1, 1, Count]), repmat(R, [1, 1, Count]));

% %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[AA_table, cc_table] = SRD_CloseLoop_GenerateTable(A_table, B_table, c_table, K_table, x_table, u_table);

% %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ode_fnc_handle = SRD_get_ode_fnc_from_ClosedLoopLinearSystem(AA_table, cc_table, time_table);


ode_fnc_handle = SRD_get_ode_fnc_from__dynamics_Linearized_Model(...
    'Handler_dynamics_generalized_coordinates_model', Handler_dynamics_generalized_coordinates_model, ...
    'Handler_dynamics_Linearized_Model', Handler_dynamics_Linearized_Model, ...
    'K_table', K_table, ...
    'x_table', x_table, ...
    'u_table', u_table, ...
    'time_table', time_table);



x0 = x_table(:, 1);
[time_table_0, solution_tape] = ode45(ode_fnc_handle, [0, tf], x0);

% figure('Color', 'w')
% plot(time_table_0, solution_tape, 'LineWidth', 3); hold on;
% plot(time_table', x_table', '--', 'LineWidth', 1);
% 
% % figure('Color', 'w')
% % plot(time_table, cc_table, 'LineWidth', 3); hold on;
% 
% 
% 
% figure('Color', 'w')

Chain = SRD_get('Chain');
DrawRobot_function = SRD_DrawRobot_get_function('DrawRobot_Type', 'STL', ... %'Default' or 'STL' or 'Custom'
    'DrawRobot_Custom_handle', [], ...
    'Function_Type', 'DrawGivenPosition', ... %'DrawGivenPosition' or 'DrawInitialPosition'  or 'DrawCurrentPosition'
    'Chain', Chain);

% DrawRobot_function(InitialPosition, [])
SRD__make_default_scene('STL')

SRD__animate__vanilla('TimeLog', time_table, 'PositionLog', x_table(1:n, :)', 'AnimationTimeLog', 0:0.05:tf, ...
    'Type', 'STL', 'DrawRobot_function', DrawRobot_function);











