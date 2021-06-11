close all; clear classes; clc;

tol = 10^(-5);

% [0, 0, phi_1, phi_2]
q = [0; 0; 0; 0];
v = zeros(size(q));
x = [q; v];
u = 0;

Handler_dynamics_generalized_coordinates_model = SRD_get('Handler_dynamics_generalized_coordinates_model');
Handler_dynamics_Linearized_Model              = SRD_get('Handler_dynamics_Linearized_Model');
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
g = f0 - A * x - B * u;


C_case0 = [...
    % 1 0 0 0 0 0 0 0;
    % 0 1 0 0 0 0 0 0;
    % 0 0 1 0 0 0 0 0;
     0 0 0 1 0 0 0 0;
    % 0 0 0 0 1 0 0 0;
    % 0 0 0 0 0 1 0 0;
    % 0 0 0 0 0 0 1 0;
    % 0 0 0 0 0 0 0 1...
    ]; 
 
C_case1 = [...
     %1 0 0 0 0 0 0 0;
     %0 1 0 0 0 0 0 0;
     0 0 1 0 0 0 0 0;
     0 0 0 1 0 0 0 0;
    % 0 0 0 0 1 0 0 0;
    % 0 0 0 0 0 1 0 0;
     0 0 0 0 0 0 1 0;
     0 0 0 0 0 0 0 1]; 
 
C_case2 = [...
     1 0 0 0 0 0 0 0;
     0 1 0 0 0 0 0 0;
     0 0 1 0 0 0 0 0;
     0 0 0 1 0 0 0 0;
    % 0 0 0 0 1 0 0 0;
    % 0 0 0 0 0 1 0 0;
     0 0 0 0 0 0 1 0;
     0 0 0 0 0 0 0 1];
 
C_case3 = [...
     1 0 0 0 0 0 0 0;
     0 1 0 0 0 0 0 0;
     0 0 1 0 0 0 0 0;
     0 0 0 1 0 0 0 0;
    % 0 0 0 0 1 0 0 0;
     0 0 0 0 0 1 0 0;
    % 0 0 0 0 0 0 1 0;
     0 0 0 0 0 0 0 1]; 

C_case4 = eye(2*n); 
 

F  = Handler_Constraints_Model.get_Jacobian(q);
dF = Handler_Constraints_Model.get_Jacobian_derivative(q);
G = [zeros(k, n), F; F, dF];
N = null(G); R = orth(G'); %E = [N, R];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%D1*dx/dt + D2*x = 0
D1 = [eye(n), zeros(n)];
D2 = [zeros(n), eye(n)];

GG = [G, zeros(2*k, 2*n); D1, D2];
GG = svd_suit(GG, tol);

V = GG.null; V = V((2*n+1):end, :);

Ra = orth(V*pinv(V)*R);
%Rn = orth((eye(2*n) - V*pinv(V))*R);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NA = svd_suit(N'*A);
Re = NA.row_space;

temp = svd_suit( (Re*Re' * V*pinv(V) * R), tol);
Rae = temp.orth;

temp2 = svd_suit( ( pinv(N'*A)*(N'*A) * V*pinv(V) * R), tol);
Rae2 = temp2.orth;


% our_R = Rae; C = C_case0; %works
% our_R = Rae; C = C_case1; %works
% our_R = Rae; C = C_case2; %works
% our_R = Rae; C = C_case3; %works
% our_R = Rae; C = C_case4; %works

% our_R = Ra; C = C_case0; % doesn't work - observer is not stabilizable
% our_R = Ra; C = C_case1; % doesn't work - observer is not stabilizable
% our_R = Ra; C = C_case2; % works
% our_R = Ra; C = C_case3; % works
% our_R = Ra; C = C_case4; % works

% our_R = R; C = C_case0; % doesn't work - observer is not stabilizable
% our_R = R; C = C_case1; % doesn't work - observer is not stabilizable
% our_R = R; C = C_case2; % doesn't work - observer is not stabilizable
our_R = R; C = C_case3; % works
% our_R = R; C = C_case4; % works






nn_1 = size(N,  2);
nn_2 = size(our_R, 2);
nn = nn_1+nn_2;

E = [N, our_R];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mm = size(C, 1);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% u_des and affine term

if norm( (eye(nn_1) - (N'*B)*pinv(N'*B)) *N'*g ) < tol
    u_des = -pinv(N'*B) *N'*g;
else
    warning('cannot create a node');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Kz    = lqr(N'*A*N, N'*B, 1*eye(nn_1), 1*eye(m));
Kzeta = pinv(N'*B) * N'*A*our_R;
K = [Kz, Kzeta];

A1 = [N'*A*N, N'*A*our_R; zeros(nn_2, nn_1), zeros(nn_2, nn_2)];
B1 = [N'*B; zeros(nn_2, m)];

L = lqr(A1', E'*C', eye(nn), eye(mm));
L = L';

CL_zLTI = [N'*A*N,          -N'*B*K;
           zeros(nn, nn_1), (A1 - L*C*E)];


Com1 = [N'*A*N-N'*B*Kz,   N'*B*K;
        zeros(nn, nn_1), (A1 - L*C*E)];
eig(Com1)    

Com2 = [N'*A*N, -N'*B*K;
        L*C*N,  (A1 - B1*K - L*C*E)];
eig(Com2)    



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PLTI in z-xi coordinates

CL_PLTI = [N'*A*N,   -N'*B*K;
           L*C*N,     (A1 - B1*K - L*C*E)];
      
CL_PLTI_g = @(zeta) [N'*A*our_R*zeta + N'*B*u_des + N'*g;
                     L*C*our_R*zeta];
                 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% LTI in x-xi coordinates
FF = R;
M = [eye(2*n),       zeros(2*n, nn),   -FF;
     zeros(nn, 2*n), eye(nn, nn),       zeros(nn, 2*k);
     G,              zeros(2*k, nn),    zeros(2*k, 2*k)];
iM = pinv(M);
iM11 = iM(1:(2*n+nn), 1:(2*n+nn));

CL_LTI_0 = [A,   -B*K;
            L*C, (A1 - B1*K - L*C*E)];
      
CL_LTI = iM11*CL_LTI_0; 
% [sort(real(eig(CL_LTI_v1)), 'descend'), sort(real(eig(CL_LTI_v2)), 'descend')]

CL_LTI_g0 = [B*u_des + g;
            zeros(nn, 1)];
CL_LTI_g = iM11*CL_LTI_g0; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% x = 0.1*randn(2*n, 1); z = N'*x; zeta = R'*x; xi_est = 0.1*randn(2*n, 1); Y0 = [z; xi_est]; 
radius = 0.01;
x = radius*[0;0;randn;randn; 0;0;0;0]; z = N'*x; zeta = our_R'*x; xi_est = radius*randn(nn, 1); Y0 = [z; xi_est]; 

g = [N'*A*our_R; L*C*our_R] * zeta;
ODEFUN = get_fnc(CL_PLTI, CL_PLTI_g(zeta));

[TOUT,YOUT] = ode45(ODEFUN,[0 25], Y0);

z        = YOUT(:, 1:nn_1);
z_est    = YOUT(:, (1+nn_1):(2*nn_1));
zeta_est = YOUT(:, (1+2*nn_1):(2*nn_1+nn_2) );
x_est = (N*z_est' + our_R*zeta_est')';

figure('Color', 'w')
subplot(2, 2, 1)
plot(TOUT, z); title('z')
subplot(2, 2, 2)
plot(TOUT, z_est); title('z est')
subplot(2, 2, 3)
if ~isempty(zeta_est)
    plot(TOUT,zeta_est); title('zeta est')
end

subplot(2, 2, 4)
plot(TOUT, x_est); title('x est')

drawnow;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ODEFUN_CLTI = get_fnc_dynamics(CL_LTI, CL_LTI_g);
Y0 = [x; xi_est]; 

[TOUT, YOUT] = ode45(ODEFUN_CLTI,[0 25], Y0);

x        = YOUT(:, 1:(2*n));
z_est    = YOUT(:, (1+2*n):(2*n+nn_1));
zeta_est = YOUT(:, (1+2*n+nn_1):(2*n+nn_1+nn_2) );
z_calc = (N'*x')';

figure('Color', 'w')
subplot(3, 2, 1)
plot(TOUT, x); title('x')
subplot(3, 2, 2)
plot(TOUT, z_est); title('z est')
subplot(3, 2, 3)

if ~isempty(zeta_est)
    plot(TOUT, zeta_est); title('zeta est')
end

subplot(3, 2, 4)
plot(TOUT, z_calc); title('z calc')
subplot(3, 2, 5)
plot(TOUT, z-z_calc); title('z-z calc')

drawnow;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function fnc = get_fnc(CL_PLTI, CL_PLTI_g)

    fnc = @my_fnc;
    function dy = my_fnc(~, y)
        dy = CL_PLTI * y + CL_PLTI_g;
    end

end

function fnc = get_fnc_dynamics(CL_LTI, CL_LTI_g)

    fnc = @my_fnc;
    function dy = my_fnc(~, y)
        dy = CL_LTI * y + CL_LTI_g;
    end

end

