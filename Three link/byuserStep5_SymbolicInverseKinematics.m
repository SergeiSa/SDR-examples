close all; clear; %clear classes;
clc; 

%Create user interfase object for SRD
SRD = SRDuserinterface;
SymbolicEngine = SRD.GetSymbolicEngine();

%%%%%%%%%%%%
%construct inverse kinematics task
rC = SymbolicEngine.GetCoM;

abs_phi = SymbolicEngine.q(1) + SymbolicEngine.q(2) + SymbolicEngine.q(3);
  
Task = [rC(1); rC(3); abs_phi]; 
%%%%%%%%%

%Call function for creating SRD inverse kinematics engine
SRD.SetupSymbolicInverseKinematics(Task);