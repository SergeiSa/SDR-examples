function c = g_dynamics_c(in1,in2)
%G_DYNAMICS_C
%    C = G_DYNAMICS_C(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 8.2.
%    15-Apr-2021 16:16:51

q1 = in1(1,:);
q2 = in1(2,:);
v1 = in2(1,:);
v2 = in2(2,:);
t2 = q1+q2;
t3 = sin(t2);
t4 = t3.*(4.9e1./2.0);
t5 = sin(q2);
c = [t4-v1+sin(q1).*(1.47e2./2.0)-t5.*v2.*(v1.*2.0+v2).*(5.0./4.0);t4-v2+t5.*v1.^2.*(5.0./4.0)];
