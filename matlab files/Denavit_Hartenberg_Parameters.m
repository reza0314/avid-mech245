clc
clear all

syms t
syms theta1(t) theta2(t) theta3(t) theta4(t) theta5(t) theta6(t)
syms dth1 dth2 dth3 dth4 dth5 dth6
syms d1 d3 d4 d6 a1 a2
syms F1 F2 F3

%% Robot parameters
% d1 = ;
% d3 = ;
% d4 = ;
% d6 = ;
%
% a1 = ;
% a2 = ;

%% D-H parameters
theta = [theta1(t) theta2(t) theta3(t) theta4(t) theta5(t) theta6(t)];
alpha = [pi/2 0 pi/2 -pi/2 pi/2 0];
d = [d1 0 -d3 d4 0 d6];
a = [a1 a2 0 0 0 0];

%% Calculate transformation matrix
T01_c = [cos(theta(1)) -cos(alpha(1))*sin(theta(1)) sin(alpha(1))*sin(theta(1)) a(1)*cos(theta(1))
    sin(theta(1)) cos(alpha(1))*cos(theta(1)) -sin(alpha(1))*cos(theta(1)) a(1)*sin(theta(1))
    0 sin(alpha(1)) cos(alpha(1)) d(1)
    0 0 0 1];
T12_c = [cos(theta(2)) -cos(alpha(2))*sin(theta(2)) sin(alpha(2))*sin(theta(2)) a(2)*cos(theta(2))
    sin(theta(2)) cos(alpha(2))*cos(theta(2)) -sin(alpha(2))*cos(theta(2)) a(2)*sin(theta(2))
    0 sin(alpha(2)) cos(alpha(2)) d(2)
    0 0 0 1];
T23_c = [cos(theta(3)) -cos(alpha(3))*sin(theta(3)) sin(alpha(3))*sin(theta(3)) a(3)*cos(theta(3))
    sin(theta(3)) cos(alpha(3))*cos(theta(3)) -sin(alpha(3))*cos(theta(3)) a(3)*sin(theta(3))
    0 sin(alpha(3)) cos(alpha(3)) d(3)
    0 0 0 1];
T34_c = [cos(theta(4)) -cos(alpha(4))*sin(theta(4)) sin(alpha(4))*sin(theta(4)) a(4)*cos(theta(4))
    sin(theta(4)) cos(alpha(4))*cos(theta(4)) -sin(alpha(4))*cos(theta(4)) a(4)*sin(theta(4))
    0 sin(alpha(4)) cos(alpha(4)) d(4)
    0 0 0 1];
T45_c = [cos(theta(5)) -cos(alpha(5))*sin(theta(5)) sin(alpha(5))*sin(theta(5)) a(5)*cos(theta(5))
    sin(theta(5)) cos(alpha(5))*cos(theta(5)) -sin(alpha(5))*cos(theta(5)) a(5)*sin(theta(5))
    0 sin(alpha(5)) cos(alpha(5)) d(5)
    0 0 0 1];
T56_c = [cos(theta(6)) -cos(alpha(6))*sin(theta(6)) sin(alpha(6))*sin(theta(6)) a(6)*cos(theta(6))
    sin(theta(6)) cos(alpha(6))*cos(theta(6)) -sin(alpha(6))*cos(theta(6)) a(6)*sin(theta(6))
    0 sin(alpha(6)) cos(alpha(6)) d(6)
    0 0 0 1];
T06_c = T01_c*T12_c*T23_c*T34_c*T45_c*T56_c;

%% Transformation matrix
T01 = [cos(theta1) 0 sin(theta1) a1*cos(theta1)
    sin(theta1) 0 -cos(theta1) a1*sin(theta1)
    0 1 0 d1
    0 0 0 1];
T12 = [cos(theta2) -sin(theta2) 0 a2*cos(theta2)
    sin(theta2) cos(theta2) 0 a2*sin(theta2)
    0 0 1 0
    0 0 0 1];
T23 = [cos(theta3) 0 sin(theta3) 0
    sin(theta3) 0 -cos(theta3) 0
    0 1 0 -d3
    0 0 0 1];
T34 = [cos(theta4) 0 -sin(theta4) 0
    sin(theta4) 0 cos(theta4) 0
    0 -1 0 d4
    0 0 0 1];
T45 = [cos(theta5) 0 sin(theta5) 0
    sin(theta5) 0 -cos(theta5) 0
    0 1 0 0
    0 0 0 1];
T56 = [cos(theta6) -sin(theta6) 0 0
    sin(theta6) cos(theta6) 0 0
    0 0 1 d6
    0 0 0 1];
T06 = T01*T12*T23*T34*T45*T56;
T06 = combine(T06(t),'sincos');

%% n
nx = T06(1,1);
ny = T06(2,1);
nz = T06(3,1);
n = [nx; ny; nz];


%% s
sx = T06(1,2);
sy = T06(2,2);
sz = T06(3,2);
s = [sx; sy; sz];

%% a
ax = T06(1,3);
ay = T06(2,3);
az = T06(3,3);
a = [ax; ay; az];

%% p
px = T06(1,4);
py = T06(2,4);
pz = T06(3,4);
p = [px; py; pz];

disp('px =')
disp(px)
disp('or')
pretty(px)

disp('py =')
disp(py)
disp('or')
pretty(py)

disp('pz =')
disp(pz)
disp('or')
pretty(pz)

%% p dot
dpx = diff(px,t);
dpy = diff(py,t);
dpz = diff(pz,t);

disp('px dot =')
disp(dpx)

disp('py dot =')
disp(dpy)

disp('pz dot =')
disp(dpz)

%% J
dpx_new = subs(dpx,[diff(theta1(t), t) diff(theta2(t), t) diff(theta3(t), t) diff(theta4(t), t) diff(theta5(t), t) diff(theta6(t), t)],[dth1 dth2 dth3 dth4 dth5 dth6]);
dpy_new = subs(dpy,[diff(theta1(t), t) diff(theta2(t), t) diff(theta3(t), t) diff(theta4(t), t) diff(theta5(t), t) diff(theta6(t), t)],[dth1 dth2 dth3 dth4 dth5 dth6]);
dpz_new = subs(dpz,[diff(theta1(t), t) diff(theta2(t), t) diff(theta3(t), t) diff(theta4(t), t) diff(theta5(t), t) diff(theta6(t), t)],[dth1 dth2 dth3 dth4 dth5 dth6]);
eqns = [dpx_new==0, dpy_new==0, dpz_new==0];
vars = [dth1 dth2 dth3 dth4 dth5 dth6];
J = equationsToMatrix(eqns,vars);
disp('J =')
disp(J)

%% Torques of joints
F = [F1; F2; F3];
tau = J.'*F;
disp('tau 1 =')
disp(tau(1))

disp('tau 2 =')
disp(tau(2))

disp('tau 3 =')
disp(tau(3))

disp('tau 4 =')
disp(tau(4))

disp('tau 5 =')
disp(tau(5))

disp('tau 6 =')
disp(tau(6))