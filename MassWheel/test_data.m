close all; clc; clear;
T = readtable('data_static.csv');

% X_ar[0,i] = state1.time
% X_ar[1:5,i] = X.reshape((1,4))
% X_ar[6,i] = state1.current
% X_ar[7,i] = state2.current
% X_ar[8,i] = u_fb
% X_ar[9:11,i] = imu_data
% X_ar[11:15,i] = X_obs[:,0]


time     = T.Var1;
x        = [T.Var2, T.Var3, T.Var4, T.Var5];
current  = [T.Var6, T.Var7];
u        = T.Var7;
u_fb     = T.Var8;
imu_data = [T.Var9, T.Var10];
x_obs    = [T.Var11, T.Var12, T.Var13, T.Var14];
x_Oleg   = [T.Var9, T.Var3, T.Var10, T.Var5];

phi2 = T.Var3;
d_phi2 = T.Var5;

figure;
subplot(2, 2, 1)
plot(time, x); hold on; title('x');
subplot(2, 2, 2)
plot(time, x_obs); hold on; title('x obs');
subplot(2, 2, 3)
plot(time, x_Oleg); hold on; title('x Oleg');


figure;
subplot(2, 2, 1)
plot(time, phi2); hold on; title('$$\phi_2$$', 'interpreter', 'latex');
subplot(2, 2, 2)
plot(time, d_phi2); hold on; title('$$\dot \phi_2$$', 'interpreter', 'latex');


figure;
subplot(2, 2, 1)
plot(1:length(time), time); hold on; title('time');
subplot(2, 2, 2)
plot(1:length(diff(time)), diff(time)); hold on; title('dt');