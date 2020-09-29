close all; clear; %clear classes; %clc;

%Create user interfase object for SRD
SRD = SRDuserinterface;
%Create ground link
Ground = SRD.CreateGroundLink();

%Create other links
index = 1;
Shin =  SRDLinkWithJoint('JointType', 'abs_pivotY', ...
                         'Order', index, ...
                         'FileName', 'datafile_Shin', ...
                         'LinkParametersStructure', [],  ...
                         'ParentLink', Ground, ...
                         'ParentFollowerNumber', 1); index = index + 1;
                     
Hip =  SRDLinkWithJoint('JointType', 'abs_pivotY', ...
                         'Order', index, ...
                         'FileName', 'datafile_Hip', ...
                         'LinkParametersStructure', [],  ...
                         'ParentLink', Shin, ...
                         'ParentFollowerNumber', 1); index = index + 1;
                     
Torso =  SRDLinkWithJoint('JointType', 'abs_pivotY', ...
                         'Order', index, ...
                         'FileName', 'datafile_Torso', ...
                         'LinkParametersStructure', [],  ...
                         'ParentLink', Hip, ...
                         'ParentFollowerNumber', 1); index = index + 1;

%Set generalised coordinates; Use negative numbers to make clock-wise
%rotations be positive
Shin.SetUsedGenCoordinates(1);
Hip.SetUsedGenCoordinates(2);
Torso.SetUsedGenCoordinates(3);

LinkArray = [Ground; Shin; Hip; Torso]; %Create array of links
InitialPosition = [pi/4; -2*pi/3; 1*pi/5]; %Define initial position of the robot

SRD.CreateRobotStructure('LinkArray', LinkArray,  ...
    'InitialPosition', InitialPosition); %Create the robot

SRD.SetupVisuals('AxisLimits', [-1; 1; -1; 1; -1; 1], ...
    'ViewAngle', [-37.5, 30], ...
    'ToDrawMeshes', false, ...
    'DrawRobot_Default_RobotColor', [0.6 0.3 0], ...
    'DrawRobot_Default_EdgeAlpha', 0.3, ...
    'DrawRobot_Default_FaceAlpha', 1, ...
    'DrawRobot_Default_LineWidth', 0.5, ...
    'DrawRobot_STL_FaceColor', [0.8 0.8 1.0], ...
    'DrawRobot_STL_EdgeColor', 'none', ...
    'DrawRobot_STL_FaceLighting', 'gouraud', ...
    'DrawRobot_STL_AmbientStrength', 0.15, ...
    'ToDrawFrames', false, ...
    'DrawRobot_Frame_Scale', 0.2, ...
    'DrawRobot_Frame_LineWidth', 1, ...
    'FileName_visuals_config', []);

DrawRobot_function = SRD_DrawRobot_get_function('DrawRobot_Type', 'Default', ... %'Default' or 'STL' or 'Custom'
    'DrawRobot_Custom_handle', [], ...
    'Function_Type', 'DrawInitialPosition', ... %'DrawGivenPosition' or 'DrawInitialPosition'  or 'DrawCurrentPosition'
    'SimulationEngine', [], ...
    'FileName_visuals_config', []);

DrawRobot_function([])


