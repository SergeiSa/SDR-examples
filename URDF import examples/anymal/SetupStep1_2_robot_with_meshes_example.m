close all; clear;
SRD = SRDuserinterface;

% LinkArray = UP_GetLinkArrayFromURDF('UrdfFilePath','./laika/urdf/laikago.urdf','ParseSTL',true);
% InitialPosition = zeros(1,21);
% 
% LinkArray = UP_GetLinkArrayFromURDF('UrdfFilePath','./iiwa/iiwa14.urdf','ParseSTL',true);
% InitialPosition = [ pi/2; pi/2; pi/2; pi/2; pi/2; pi/2; pi/2; pi/2; pi/2; 0; 0];
% 
% AxisLimits = []; %Set axis limits for the snapshot
% ViewAngle = [45, 45];
% DrawFrames = true;
% DrawMeshes = true;
% SRD.CreateRobotStructure(LinkArray, InitialPosition, AxisLimits, ViewAngle,DrawFrames,DrawMeshes); %Create the robot
% axis equal;

%Set generalised coordinates; Use negative numbers to make clock-wise
%rotations be positive
% Shin.SetUsedGenCoordinates(1);
% Hip.SetUsedGenCoordinates(2);
% Torso.SetUsedGenCoordinates(3);


LinkArray = SRD_GetLinkArrayFromURDF('UrdfFilePath','./urdf/anymal.urdf','ParseSTL',true);

Color_0 = [1, 0, 0];
Color_1 = [0.95, 0.5, 0.3];
Color_2 = [0.00, 0.00, 0.00];
SRD_LinkSet_Color(LinkArray, 'base', Color_0);
SRD_LinkSet_Color(LinkArray, 'base_inertia', Color_0);
SRD_LinkSet_Color(LinkArray, 'LF_HIP', Color_1);
SRD_LinkSet_Color(LinkArray, 'LF_THIGH', Color_1);
SRD_LinkSet_Color(LinkArray, 'LF_SHANK', Color_2);
SRD_LinkSet_Color(LinkArray, 'LF_ADAPTER', Color_1);
SRD_LinkSet_Color(LinkArray, 'LF_FOOT', Color_2);
SRD_LinkSet_Color(LinkArray, 'RF_HIP', Color_1);
SRD_LinkSet_Color(LinkArray, 'RF_THIGH', Color_1);
SRD_LinkSet_Color(LinkArray, 'RF_SHANK', Color_2);
SRD_LinkSet_Color(LinkArray, 'RF_ADAPTER', Color_1);
SRD_LinkSet_Color(LinkArray, 'RF_FOOT', Color_2);
SRD_LinkSet_Color(LinkArray, 'LH_HIP', Color_1);
SRD_LinkSet_Color(LinkArray, 'LH_THIGH', Color_1);
SRD_LinkSet_Color(LinkArray, 'LH_SHANK', Color_2);
SRD_LinkSet_Color(LinkArray, 'LH_ADAPTER', Color_1);
SRD_LinkSet_Color(LinkArray, 'LH_FOOT', Color_2);
SRD_LinkSet_Color(LinkArray, 'RH_HIP', Color_1);
SRD_LinkSet_Color(LinkArray, 'RH_THIGH', Color_1);
SRD_LinkSet_Color(LinkArray, 'RH_SHANK', Color_2);
SRD_LinkSet_Color(LinkArray, 'RH_ADAPTER', Color_1);
SRD_LinkSet_Color(LinkArray, 'RH_FOOT', Color_2);

% InitialPosition = [0, 0, pi/2,0, 0,pi/2, 0,0,pi/2,0,0,pi/2]';
InitialPosition = zeros(12,1);
SRD_save(LinkArray, 'LinkArray');
SRD_save(InitialPosition, 'InitialPosition');

Chain = SRD_Chain(LinkArray);
SRD_save(Chain, 'Chain');
%Chain.Update(InitialPosition)


SRD_SetupVisuals('AxisLimits', [-1; 1; -1; 1; -1; 1], ...
    'ViewAngle', [-37.5, 30], ...
    'ToDrawMeshes', true, ...
    'Animation_ToUseGrid', true, ...
    'Animation_ToUseGridMinor', true, ...
    'DrawRobot_Default_RobotColor', [0.6 0.3 0], ...
    'DrawRobot_Default_EdgeAlpha', 0.3, ...
    'DrawRobot_Default_FaceAlpha', 1, ...
    'DrawRobot_Default_LineWidth', 0.5, ...
    'DrawRobot_STL_FaceColor', [0.8 0.8 1.0], ...
    'DrawRobot_STL_EdgeColor', 'none', ...
    'DrawRobot_STL_FaceLighting', 'gouraud', ...
    'DrawRobot_STL_AmbientStrength', 0.15, ...
    'DrawRobot_STL_camlight', 'headlight', ...
    'DrawRobot_STL_material', 'metal', ... %shiny dull metal
    'ToDrawFrames', false, ...
    'DrawRobot_Frame_Scale', 0.2, ...
    'DrawRobot_Frame_LineWidth', 1);

DrawRobot_function = SRD_DrawRobot_get_function('DrawRobot_Type', 'STL', ... %'Default' or 'STL' or 'Custom'
    'DrawRobot_Custom_handle', [], ...
    'Function_Type', 'DrawGivenPosition', ... %'DrawGivenPosition' or 'DrawInitialPosition'  or 'DrawCurrentPosition'
    'Chain', Chain);

DrawRobot_function(InitialPosition, [])
SRD__make_default_scene('STL')
            