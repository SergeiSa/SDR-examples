function Task = g_Constraint(in1)
%G_CONSTRAINT
%    TASK = G_CONSTRAINT(IN1)

%    This function was generated by the Symbolic Math Toolbox version 8.2.
%    31-Mar-2021 12:05:44

q1 = in1(1,:);
q2 = in1(2,:);
q3 = in1(3,:);
t2 = q1+q2+q3;
t3 = cos(t2);
t4 = q1+q2;
t5 = cos(t4);
t6 = cos(q1);
Task = [t3./2.0+t5./2.0+t6./2.0;t3+t5+t6];
