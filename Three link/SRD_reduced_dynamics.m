function description = SRD_reduced_dynamics(varargin)
Parser = inputParser;
Parser.FunctionName = 'SRD_reduced_dynamics';
% Parser.addOptional('SymbolicEngine', []);

Parser.addOptional('Symbolic_ToOptimizeFunctions', true);
Parser.addOptional('Casadi_cfile_name', 'transv_dynamics');

%H*ddq + C*dq + g = T*u
Parser.addOptional('H_callable', []);
Parser.addOptional('C_callable', []);
Parser.addOptional('g_callable', []);
Parser.addOptional('T_callable', []);


Parser.addOptional('FunctionName_alpha', 'get_alpha');
Parser.addOptional('FunctionName_beta', 'get_beta');
Parser.addOptional('FunctionName_gamma', 'get_gamma');

Parser.addOptional('Path', 'autogen_reduced_dynamics/');

Parser.addOptional('N_dof', []);

Parser.addOptional('c0', []); % s = c0*q
Parser.addOptional('H0', []); % H0*q = Phi(s)

Parser.parse(varargin{:});



% if isempty(Parser.Results.SymbolicEngine)
%     error('Please provide SymbolicEngine')
% else
%     SymbolicEngine = Parser.Results.SymbolicEngine;
% end

if ~isempty(Parser.Results.Path)
    if ~exist(Parser.Results.Path, 'dir')
        mkdir(Parser.Results.Path)
    end
end

N_dof = Parser.Results.N_dof;
c0 = Parser.Results.c0;
H0 = Parser.Results.H0;

% Generating function for reduced dynamics of the pednubot

% -----------------------------------------------------------------------
% Symbolic variables for virtual constraints
% -----------------------------------------------------------------------
Phi = sym('phi', [N_dof,1], 'real');
Phi_prm = sym('phi%d_prm', [N_dof,1], 'real');
Phi_2prm = sym('phi%d_2prm', [N_dof,1], 'real');

% ----------------------------------------------------------------------
% Compute reduced dynamics
% ----------------------------------------------------------------------
% Annihilator of B matrix
B = Parser.Results.T_callable();
B_anh = null(B')';

% Change of coordinates from q to Phi
M_phi = Parser.Results.H_callable(Phi);
C_phi = Parser.Results.C_callable(Phi, Phi_prm); 
G_phi = Parser.Results.g_callable(Phi);

% Obtaining alpha, beta, gamma expressions
alpha = B_anh * M_phi * Phi_prm;
beta = B_anh * (M_phi * Phi_2prm + C_phi * Phi_prm );
gamma = B_anh * G_phi;


% Generate functions
% -------------------------------------------------------------------
FileName_alpha = [Parser.Results.Path, Parser.Results.FunctionName_alpha];
FileName_beta = [Parser.Results.Path, Parser.Results.FunctionName_beta];
FileName_gamma = [Parser.Results.Path, Parser.Results.FunctionName_gamma];
        
        
matlabFunction(alpha, 'File', FileName_alpha, 'Vars', {Phi, Phi_prm});
matlabFunction(beta, 'File', FileName_beta, 'Vars', {Phi, Phi_prm, Phi_2prm});
matlabFunction(gamma, 'File', FileName_gamma, 'Vars', {Phi});


description.c0 = c0;  
description.H0 = H0;   


description.FileName_alpha = FileName_alpha;   
description.FileName_beta  = FileName_beta;  
description.FileName_gamma = FileName_gamma;  


end