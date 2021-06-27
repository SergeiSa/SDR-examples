close all; clear classes; clc;

InitialPosition = SRD_get('InitialPosition');

% Handler_State = SRDHandler_State(...
%     'InitialPosition', InitialPosition, ...
%     'InitialVelocity', zeros(size(InitialPosition)));

Handler_StateSpace = SRDHandler_StateSpace('InitialState', [InitialPosition; zeros(size(InitialPosition))]);

Handler_StateConverter_StateSpace2GenCoord = SRDHandler_StateConverter_StateSpace2GenCoord(...
    'Handler_StateSpace', Handler_StateSpace);
Handler_State = Handler_StateConverter_StateSpace2GenCoord;

Handler_IK_Solution = SRD_get('Handler_IK_Solution');


% %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Handler_dynamics_generalized_coordinates_model = SRD_get('Handler_dynamics_generalized_coordinates_model');
Handler_dynamics_Linearized_Model = SRD_get('Handler_dynamics_Linearized_Model');


Handler_dynamics_GC_model_evaluator = SRDHandler_dynamics_GC_model_evaluator(...
    'Handler_dynamics_generalized_coordinates_model', Handler_dynamics_generalized_coordinates_model, ...
    'Handler_State', Handler_State, ...
    'UsePinv', true);

Handler_Constraints_Model = SRD_get('Handler_Constraints_Model');

Handler_dynamics_Linear_model_evaluator = SRD_get_handler__dynamics_Linear_model_evaluator(...
    'Handler_dynamics_generalized_coordinates_model', Handler_dynamics_GC_model_evaluator, ...
    'Handler_dynamics_Linearized_Model', Handler_dynamics_Linearized_Model, ...
    'Handler_State', Handler_State, ...
    'Handler_Controller', [], ...
    'ToEvaluate_c', false);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %

dt = 0.001;
tf = Handler_IK_Solution.TimeExpiration;
% tf = 0.15;

Handler_Time = SRDHandler_Time('TimeLog', 0:dt:tf);

% %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Handler_Desired_State = SRD_get_handler__desired_state_static(...
%     'static_q',  InitialPosition, ...
%     'dof_robot', Handler_dynamics_generalized_coordinates_model.dof_configuration_space_robot);

Handler_Desired_State = SRD_get_handler__desired_state(...
    'Handler_ControlInput', Handler_IK_Solution, ...
    'Handler_Time',         Handler_Time);

Handler_State_StateSpace = SRDHandler_StateConverter_GenCoord2StateSpace(...
    'Handler_State', Handler_State);

Handler_Desired_State_StateSpace = SRDHandler_StateConverter_GenCoord2StateSpace(...
    'Handler_State', Handler_Desired_State);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %


% %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Handler_InverseDynamics = SRD_get_handler__InverseDynamics_Vanilla__desired_trajectory(...
    'Handler_ControlInput', Handler_Desired_State, ...
    'Handler_dynamics_generalized_coordinates_model', Handler_dynamics_GC_model_evaluator, ...
    'Handler_Time', Handler_Time);

Handler_InverseDynamics = SRD_get_handler__InverseDynamicsConstrained_QR(...
    'Handler_ControlInput', Handler_Desired_State, ...
    'Handler_dynamics_generalized_coordinates_model', Handler_dynamics_generalized_coordinates_model, ...
    'Handler_Constraints_Model', Handler_Constraints_Model, ...
    'Handler_Time', Handler_Time);

Handler_ComputedTorqueController = SRD_get_handler__ComputedTorqueController(...
    'Handler_State', Handler_State, ...
    'Handler_ControlInput', Handler_Desired_State, ...
    'Handler_dynamics_generalized_coordinates_model', Handler_dynamics_GC_model_evaluator, ...
    'Handler_Time', Handler_Time, ...
    'Handler_InverseDynamics', Handler_InverseDynamics, ...
    'Kp', 500*eye(Handler_Desired_State.dof_robot), ...
    'Kd', 100*eye(Handler_Desired_State.dof_robot));

% Handler_LQR = SRD_get_handler__LQR_Controller(...
%     'Handler_State_StateSpace', Handler_State_StateSpace, ...
%     'Handler_ControlInput_StateSpace', Handler_Desired_State_StateSpace, ...
%     'Handler_dynamics_Linearized_Model', Handler_dynamics_Linear_model_evaluator, ...
%     'Handler_Time', Handler_Time, ...
%     'Handler_InverseDynamics', Handler_InverseDynamics, ...
%     'Q', 10*eye(Handler_dynamics_Linear_model_evaluator.dof_robot_StateSpace), ...
%     'R', 1*eye(Handler_dynamics_Linear_model_evaluator.dof_control));

Handler_LQR = SRD_get_handler__Constrained_LQR_Controller(...
    'Handler_State', Handler_State, ...
    'Handler_State_StateSpace', Handler_State_StateSpace, ...
    'Handler_Constraints_Model', Handler_Constraints_Model, ...
    'Handler_ControlInput_StateSpace', Handler_Desired_State_StateSpace, ...
    'Handler_dynamics_Linearized_Model', Handler_dynamics_Linear_model_evaluator, ...
    'Handler_Time', Handler_Time, ...
    'Handler_InverseDynamics', Handler_InverseDynamics, ...
    'Q', 10*eye(Handler_dynamics_Linear_model_evaluator.dof_robot_StateSpace), ...
    'R', 1*eye(Handler_dynamics_Linear_model_evaluator.dof_control));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %

% MainController = Handler_ComputedTorqueController;
MainController = Handler_LQR;

Handler_dynamics_Linear_model_evaluator.Handler_Controller = Handler_InverseDynamics;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %


% Handler_solver_Taylor = SRD_get_handler__solver_Taylor(...
%     'Handler_State', Handler_State, ...
%     'Handler_Controller', MainController, ...
%     'Handler_dynamics_generalized_coordinates_model', Handler_dynamics_GC_model_evaluator, ...
%     'Handler_Time', Handler_Time);

Handler_solver_Euler = SRD_get_handler__solver_Euler(...
    'Handler_StateSpace', Handler_StateSpace, ...
    'Handler_Controller', MainController, ...
    'Handler_FirstOrderSystem', @(x, u) Handler_dynamics_GC_model_evaluator.get_FirstOrderSystem_qv(x, u), ...
    'Handler_Time', Handler_Time);


Handler_solver_TaylorConstrained = SRD_get_handler__solver_TaylorConstrained(...
    'Handler_State', Handler_State, ...
    'Handler_Controller', MainController, ...
    'Handler_dynamics_generalized_coordinates_model', Handler_dynamics_GC_model_evaluator, ...
    'Handler_Time', Handler_Time, ...
    'Handler_Constraints_Model', Handler_Constraints_Model);



% %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Handler_State_Logger_vanilla = SRD_get_handler__State_Logger__vanilla(...
    'Handler_State', Handler_State, ...
    'Handler_Time', Handler_Time, ...
    'ToLogAcceleration',  false);

Handler_SimulationTickDisplay = SRD_get_handler__SimulationTickDisplay(...
    'Handler_Simulation', Handler_Time, ...
    'DisplayOneTickIn', 100);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %

Handler_Updater = SRDHandler_Updater({...
    Handler_StateConverter_StateSpace2GenCoord, ...
    Handler_Desired_State, ...
    Handler_State_StateSpace, ...
    Handler_Desired_State_StateSpace, ...
    Handler_dynamics_GC_model_evaluator,...
    Handler_InverseDynamics, ...
    Handler_dynamics_Linear_model_evaluator, ...
    Handler_LQR, ...
    Handler_solver_TaylorConstrained, ...
    Handler_State_Logger_vanilla, ...
    Handler_SimulationTickDisplay, ...
    Handler_Time ...
    });

% Handler_Desired_State, Handler_State_StateSpace, Handler_Desired_State_StateSpace, ...
%     Handler_dynamics_GC_model_evaluator;
% Handler_InverseDynamics, Handler_ComputedTorqueController;
% Handler_solver_Taylor;


for i = 1:(length(Handler_Time.TimeLog)-1)
    Handler_Updater.Update();
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


figure('Color', 'w', 'Name', 'Positions');

subplot(1, 2, 1)
SRDgraphic_PlotGeneric(Handler_Time.TimeLog, Handler_State_Logger_vanilla.Log.q, ...
    'NewFigure', false, 'FigureName', 'Generic', ...
    'LableVariable', 'q', 'Title', []);

subplot(1, 2, 2)
SRDgraphic_PlotGeneric(Handler_Time.TimeLog, Handler_State_Logger_vanilla.Log.v, ...
    'NewFigure', false, 'FigureName', 'Generic', ...
    'LableVariable', 'v', 'Title', []);

drawnow;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

ToAnimate = false;

if ToAnimate
DrawRobot_function = SRD_DrawRobot_get_function('DrawRobot_Type', 'Default', ... %'Default' or 'STL' or 'Custom'
    'DrawRobot_Custom_handle', [], ...
    'Function_Type', 'DrawGivenPosition', ... %'DrawGivenPosition' or 'DrawInitialPosition'  or 'DrawCurrentPosition'
    'SimulationEngine', [], ...
    'FileName_visuals_config', []); %use default visuals

SRD__animate__vanilla('Handler_Simulation', Handler_Simulation, ...
    'Handler_Logger', Handler_State_Logger_vanilla, ...
    'AnimationTimeLog', 0:10*dt:tf, ...
    'DrawRobot_function', DrawRobot_function, ...
    'NewFigure', true, ...
    'FigureName', 'Animation', ...
    'FileName_visuals_config', []);
end



