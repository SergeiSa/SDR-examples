function Constraint = g_dynamics_LagrangeMultiplier_Constraint(in1)
%G_DYNAMICS_LAGRANGEMULTIPLIER_CONSTRAINT
%    CONSTRAINT = G_DYNAMICS_LAGRANGEMULTIPLIER_CONSTRAINT(IN1)

%    This function was generated by the Symbolic Math Toolbox version 7.0.
%    30-Oct-2017 22:52:10

q1 = in1(1,:);
q2 = in1(2,:);
q3 = in1(3,:);
q4 = in1(4,:);
q5 = in1(5,:);
q6 = in1(6,:);
q7 = in1(7,:);
q8 = in1(8,:);
q9 = in1(9,:);
q10 = in1(10,:);
q11 = in1(11,:);
t2 = sin(q4);
t3 = sin(q6);
t4 = cos(q4);
t5 = cos(q6);
t6 = sin(q5);
t7 = cos(q10);
t8 = cos(q5);
t9 = sin(q10);
t10 = sin(q7);
t11 = sin(q9);
t12 = cos(q7);
t13 = cos(q9);
t14 = sin(q8);
t15 = cos(q11);
t16 = cos(q8);
t17 = sin(q11);
Constraint = [q1-t7.*(t2.*t3+t4.*t5.*t6).*(1.0./2.0)-t2.*t3.*(1.0./2.0)-sqrt(3.0).*(1.0./4.0)-t4.*t5.*t6.*(1.0./2.0)-t5.*t8.*t9.*(1.0./2.0);q2+t7.*(t2.*t5-t3.*t4.*t6).*(1.0./2.0)+t2.*t5.*(1.0./2.0)-t3.*t4.*t6.*(1.0./2.0)-t3.*t8.*t9.*(1.0./2.0);q3-t4.*t8.*(1.0./2.0)+t6.*t9.*(1.0./2.0)-t4.*t7.*t8.*(1.0./2.0)-2.499999999999999e-1;q1-t15.*(t10.*t11+t12.*t13.*t14).*(1.0./2.0)-t10.*t11.*(1.0./2.0)-t12.*t13.*t14.*(1.0./2.0)-t13.*t16.*t17.*(1.0./2.0)+1.385729807284855e-1;q2+t15.*(t10.*t13-t11.*t12.*t14).*(1.0./2.0)+t10.*t13.*(1.0./2.0)-t11.*t12.*t14.*(1.0./2.0)-t11.*t16.*t17.*(1.0./2.0);q3-t12.*t16.*(1.0./2.0)+t14.*t17.*(1.0./2.0)-t12.*t15.*t16.*(1.0./2.0)-2.523271025828619e-1];
