function A = g_linearization_A(in1,in2,in3,in4)
%G_LINEARIZATION_A
%    A = G_LINEARIZATION_A(IN1,IN2,IN3,IN4)

%    This function was generated by the Symbolic Math Toolbox version 8.4.
%    19-Nov-2020 18:56:34

iH1_1 = in4(1);
iH1_2 = in4(4);
iH1_3 = in4(7);
iH2_1 = in4(2);
iH2_2 = in4(5);
iH2_3 = in4(8);
iH3_1 = in4(3);
iH3_2 = in4(6);
iH3_3 = in4(9);
q1 = in1(1,:);
q2 = in1(2,:);
q3 = in1(3,:);
u1 = in3(1,:);
u2 = in3(2,:);
u3 = in3(3,:);
v1 = in2(1,:);
v2 = in2(2,:);
v3 = in2(3,:);
t2 = cos(q2);
t3 = cos(q3);
t4 = sin(q1);
t5 = sin(q2);
t6 = sin(q3);
t7 = q1+q2;
t8 = q2+q3;
t9 = v1+v2;
t10 = v1.^2;
t11 = v2.^2;
t12 = v3.^2;
t19 = -u1;
t20 = -v1;
t21 = -v2;
t22 = -v3;
t13 = cos(t7);
t14 = cos(t8);
t15 = q3+t7;
t16 = sin(t7);
t17 = sin(t8);
t18 = t9+v3;
t25 = t3.*(5.0./2.0);
t26 = t6.*(5.0./2.0);
t27 = t6.*(5.0./4.0);
t28 = t5.*(1.5e+1./2.0);
t29 = t5.*(1.5e+1./4.0);
t40 = t4.*(2.45e+2./2.0);
t41 = t6.*v3.*(-5.0./2.0);
t62 = t6.*t12.*(-5.0./4.0);
t23 = cos(t15);
t24 = sin(t15);
t30 = t26.*v1;
t31 = t26.*v2;
t32 = t26.*v3;
t33 = t6+t17;
t34 = t14.*(5.0./4.0);
t35 = t17.*(5.0./2.0);
t36 = t17.*(5.0./4.0);
t42 = t28.*v1;
t43 = t28.*v2;
t45 = t10.*t27;
t46 = t11.*t27;
t47 = t12.*t27;
t48 = t41.*v1;
t49 = t41.*v2;
t51 = -t40;
t60 = t13.*(1.47e+2./2.0);
t61 = t16.*(1.47e+2./2.0);
t63 = t10.*t29;
t64 = t11.*t29;
t37 = t30.*v2;
t38 = t30.*v3;
t39 = t31.*v3;
t44 = t32+1.0;
t50 = t42.*v2;
t52 = t34.*v2;
t53 = t35.*v1;
t54 = t35.*v2;
t55 = t35.*v3;
t56 = t36.*v2;
t65 = t23.*(4.9e+1./2.0);
t66 = t24.*(4.9e+1./2.0);
t67 = -t61;
t69 = t10.*t34;
t70 = t10.*t36;
t71 = t11.*t36;
t72 = t12.*t36;
t73 = t26+t35;
t74 = t26+t36;
t75 = t27+t36;
t76 = t28+t35;
t77 = t29+t36;
t57 = t53.*v2;
t58 = t53.*v3;
t59 = t54.*v3;
t68 = -t66;
t78 = t73.*v3;
t79 = t74.*v3;
t80 = t75.*v3;
t81 = t76.*v2;
t82 = t77.*v2;
t83 = t30+t31+t53;
t84 = t41+t42+t53;
t95 = t32+t42+t43+t53+t54+t55;
t96 = t22+t37+t45+t46+t66+t70+u3;
t99 = t21+t48+t49+t61+t62+t63+t66+t70+u2;
t85 = t56+t80;
t88 = t78+t81;
t89 = t79+t82;
t97 = iH2_3.*t96;
t98 = iH3_3.*t96;
t100 = iH2_2.*t99;
t101 = iH3_2.*t99;
t102 = t19+t38+t39+t47+t50+t51+t57+t58+t59+t64+t67+t68+t71+t72+v1;
t86 = t85.*v3;
t87 = t22.*t85;
t90 = t88+1.0;
t91 = t88.*v1;
t92 = t89.*v2;
t93 = t20.*t88;
t94 = t21.*t89;
t103 = iH2_1.*t102;
t104 = iH3_1.*t102;
t105 = -t103;
t106 = -t104;
t107 = t20+t40+t61+t66+t87+t93+t94+u1;
t108 = iH3_1.*t107;
t109 = t97+t100+t105;
t110 = t98+t101+t106;
t111 = t36.*t109;
t112 = t98+t101+t108;
t113 = t73.*t110;
t114 = t74.*t110;
t115 = t75.*t110;
t116 = t76.*t109;
t117 = t77.*t109;
t118 = t111+t115;
t119 = t113+t116;
t120 = t114+t117;
A = reshape([0.0,0.0,0.0,t60+t65+cos(q1).*(2.45e+2./2.0)+iH1_1.*t119+iH1_3.*t118+iH1_2.*t120,t60+t65+iH2_1.*t119+iH2_3.*t118+iH2_2.*t120,t65+iH3_1.*t119+iH3_3.*t118+iH3_2.*t120,0.0,0.0,0.0,t60+t65+iH1_1.*t120-t2.*t11.*(1.5e+1./4.0)-t11.*t14.*(5.0./4.0)-t12.*t14.*(5.0./4.0)+iH1_2.*t26.*t110+iH1_3.*t27.*t110-t2.*v1.*v2.*(1.5e+1./2.0)-t14.*v1.*v2.*(5.0./2.0)-t14.*v1.*v3.*(5.0./2.0)-t14.*v2.*v3.*(5.0./2.0),t60+t65+t69+iH2_1.*t120+t2.*t10.*(1.5e+1./4.0)+iH2_2.*t26.*t110+iH2_3.*t27.*t110,t65+t69+iH3_1.*t120+iH3_2.*t26.*t110+iH3_3.*t27.*t110,0.0,0.0,0.0,t65+t22.*(t52+v3.*(t3.*(5.0./4.0)+t34))+t21.*(t52+v3.*(t25+t34))+iH1_1.*(t75.*t112+t36.*(t97+t100+iH2_1.*t107))+t20.*(t14.*v2.*(5.0./2.0)+v3.*(t14.*(5.0./2.0)+t25))+iH1_2.*t27.*t112,t65+t69+iH2_1.*t118-t3.*t12.*(5.0./4.0)+iH2_2.*t27.*t110-t3.*v1.*v3.*(5.0./2.0)-t3.*v2.*v3.*(5.0./2.0),t65+t69+iH3_1.*t118+t3.*t10.*(5.0./4.0)+t3.*t11.*(5.0./4.0)+iH3_2.*t27.*t110+t25.*v1.*v2,1.0,0.0,0.0,iH1_2.*t84+iH1_3.*t83-iH1_1.*t90,iH2_2.*t84+iH2_3.*t83-iH2_1.*t90,iH3_2.*t84+iH3_3.*t83-iH3_1.*t90,0.0,1.0,0.0,-iH1_2.*t44-iH1_1.*t95+iH1_3.*t9.*t26,-iH2_2.*t44-iH2_1.*t95+iH2_3.*t9.*t26,-iH3_2.*t44-iH3_1.*t95+iH3_3.*t9.*t26,0.0,0.0,1.0,-iH1_3-iH1_2.*t6.*t18.*(5.0./2.0)-iH1_1.*t18.*t33.*(5.0./2.0),-iH2_3-iH2_2.*t6.*t18.*(5.0./2.0)-iH2_1.*t18.*t33.*(5.0./2.0),-iH3_3-iH3_2.*t6.*t18.*(5.0./2.0)-iH3_1.*t18.*t33.*(5.0./2.0)],[6,6]);
