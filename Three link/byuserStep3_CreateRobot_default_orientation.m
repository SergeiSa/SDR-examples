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

%an example of how to use default orientation
Shin.PivotZeroOrientation = roty(-20);
Hip.PivotZeroOrientation = roty(-20);
Torso.PivotZeroOrientation = roty(-20);

LinkArray = [Ground; Shin; Hip; Torso]; %Create array of links
InitialPosition = [pi/4; -2*pi/3; 1*pi/5]; %Define initial position of the robot
% InitialPosition = [pi/6; -pi/6; pi/6]; %Define initial position of the robot
AxisLimits = [-1; 1; -1; 1; 0; 2]; %Set axis limits for the snapshot
ViewAngle = [0, 0];

SRD.CreateRobotStructure(LinkArray, InitialPosition, AxisLimits, ViewAngle); %Create the robot