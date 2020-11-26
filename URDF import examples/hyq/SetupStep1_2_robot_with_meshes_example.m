close all; clear;
SRD = SRDuserinterface;




LinkArray = SRD_GetLinkArrayFromURDF('UrdfFilePath','./urdf/hyq.urdf','ParseSTL',true);

Color_0 = [1, 0, 0];
Color_1 = [0.95, 0.5, 0.3];
Color_2 = [0.00, 0.00, 0.00];
Color_3 = [0.1, 0.1, 0.1];

SRD_LinkSet_Color(LinkArray, 'trunk', Color_0);
SRD_LinkSet_Color(LinkArray, 'lf_haa_joint', Color_1);
SRD_LinkSet_Color(LinkArray, 'lf_hipassembly', Color_1);
SRD_LinkSet_Color(LinkArray, 'lf_upperleg', Color_2);
SRD_LinkSet_Color(LinkArray, 'lf_lowerleg', Color_1);
SRD_LinkSet_Color(LinkArray, 'lf_foot', Color_2);
SRD_LinkSet_Color(LinkArray, 'rf_hipassembly', Color_1);
SRD_LinkSet_Color(LinkArray, 'rf_upperleg', Color_2);
SRD_LinkSet_Color(LinkArray, 'rf_lowerleg', Color_1);
SRD_LinkSet_Color(LinkArray, 'rf_foot', Color_2);
SRD_LinkSet_Color(LinkArray, 'lh_hipassembly', Color_1);
SRD_LinkSet_Color(LinkArray, 'lh_upperleg', Color_2);
SRD_LinkSet_Color(LinkArray, 'lh_lowerleg', Color_1);
SRD_LinkSet_Color(LinkArray, 'lh_foot', Color_2);
SRD_LinkSet_Color(LinkArray, 'rh_hipassembly', Color_1);
SRD_LinkSet_Color(LinkArray, 'rh_upperleg', Color_2);
SRD_LinkSet_Color(LinkArray, 'rh_lowerleg', Color_1);
SRD_LinkSet_Color(LinkArray, 'rh_foot', Color_2);
% SRD_LinkSet_Color(LinkArray, 'RL_thigh', Color_1);
% SRD_LinkSet_Color(LinkArray, 'RL_calf', Color_3);
% SRD_LinkSet_Color(LinkArray, 'RL_foot', Color_2);
% SRD_LinkSet_Color(LinkArray, 'RH_FOOT', Color_2);

% InitialPosition = [ pi/2; pi/2; pi/2; pi/2; pi/2; pi/2; pi/2; pi/2; pi/2; 0; 0];
% InitialPosition = [ pi/2; 0.2*pi/2; 0.2*pi/2; -1.0*pi/2; 1.2*pi/2; 1.2*pi/2; pi/2;];

InitialPosition = [0, 0, pi/2,0, 0,-pi/2, 0,0,pi/2,0,0,-pi/2]';

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
    'DrawRobot_Default_RobotColor', [0.6 0.3 0.1], ...
    'DrawRobot_Default_EdgeAlpha', 0.3, ...
    'DrawRobot_Default_FaceAlpha', 1, ...
    'DrawRobot_Default_LineWidth', 0.5, ...
    'DrawRobot_STL_FaceColor', [0.95 0.8 0.2], ...
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
            