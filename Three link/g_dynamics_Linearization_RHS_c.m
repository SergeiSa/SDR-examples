function out1 = g_dynamics_Linearization_RHS_c(in1,in2)
%G_DYNAMICS_LINEARIZATION_RHS_C
%    OUT1 = G_DYNAMICS_LINEARIZATION_RHS_C(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 8.1.
%    14-Nov-2018 09:10:20

s1 = in1(1,:);
s2 = in1(2,:);
s3 = in1(3,:);
out1 = [0.0;0.0;0.0;sin(s1).*2.94e2-s1.*cos(s1).*2.94e2;sin(s2).*(8.82e2./5.0)-s2.*cos(s2).*(8.82e2./5.0);sin(s3).*(2.94e2./5.0)-s3.*cos(s3).*(2.94e2./5.0)];
