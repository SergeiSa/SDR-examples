function B = g_linearization_B(in1,in2,in3)
%G_LINEARIZATION_B
%    B = G_LINEARIZATION_B(IN1,IN2,IN3)

%    This function was generated by the Symbolic Math Toolbox version 8.4.
%    19-Nov-2020 18:56:36

iH1_1 = in3(1);
iH1_2 = in3(4);
iH1_3 = in3(7);
iH2_1 = in3(2);
iH2_2 = in3(5);
iH2_3 = in3(8);
iH3_1 = in3(3);
iH3_2 = in3(6);
iH3_3 = in3(9);
B = reshape([0.0,0.0,0.0,iH1_1,iH2_1,iH3_1,0.0,0.0,0.0,iH1_2,iH2_2,iH3_2,0.0,0.0,0.0,iH1_3,iH2_3,iH3_3],[6,3]);
