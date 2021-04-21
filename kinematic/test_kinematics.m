close all;
l=0.205;
DH = [  1, 0, l, 0, 0;
        2, 0, l, 0, 0;
        3, 0, 0, 0, 0;];
DH2 =  [4, 0, l, 0, 0;
        5, 0, l, 0, 0;
        6, 0, 0, 0, 0;];
arm1 = SerialLink(DH, 'name', 'arm1');
arm1.base = transl(-0.125, 0, 0);
arm2 = SerialLink(DH2, 'name', 'arm2');
arm2.base = transl(0.125, 0, 0);

% Inverse kinematics validation
p0 = [0; 0.35];
q0 = ikine5(p0);
% show working modes
% for i=1:4
% arm1.plot(q0(1:3,i)')
% hold on;
% arm2.plot([q0(4:5,i)',0])
% pause();
% end
q0 = q0(:,4); % Using -+ mode (4)
qa0 = [q0(1);q0(4)];

% Forward kinematics validation
[p1, q1] = fkine5(qa0,1);
[p2, q2] = fkine5(qa0,-1);
%show assembly modes
% plot(p1(1),p1(2),'k*');
% hold on;
% plot(p2(1), p2(2), 'k*');
% arm1.plot(q1(1:3)')
% arm2.plot([q1(4:5)',0])
% pause();
% arm1.plot(q2(1:3)')
% arm2.plot([q2(4:5)',0])
% pause();

% first order inverse kinematics
t0r = [1;1];
dt = 0.1;
pv1 = p0 + t0r * dt;
[qd,qda,qdu,Ar,B,J,Jinv,Jt,Jta,Jtd,Psit] = ifirstkine5(t0r,p0,q0);
disp('active joint velocities');
qda
% quiver(p0(1),p0(2),pv1(1)-p0(1),pv1(2)-p0(2),0);
% hold on;
% arm1.plot(q1(1:3)')
% arm2.plot([q1(4:5)',0])
% pause();

% first order forward kinematics
[t0r2] = ffirstkine5(q0,qda);
pv2 = p0 + t0r2 * dt;
disp('platform twist');
t0r2
% quiver(p0(1),p0(2),pv2(1)-p0(1),pv2(2)-p0(2),0);
% hold on;
% arm1.plot(q1(1:3)')
% arm2.plot([q1(4:5)',0])
% pause();

% second order inverse kinematics
t0rd = [-1;1];
pa1 = p0 +   t0r * dt + t0rd * dt^2 / 2;
[qdd,qdda,qddu,Psit1,Psit2,atp,ad,b0p1,b0p2,Psidt] = isecondkine5(t0rd,t0r,p0,q0,qd,Ar,B,Jinv,Jt,Jta,Jtd,Psit);
disp('active joints acceleration');
qdda
% quiver(p0(1),p0(2),pa1(1)-p0(1),pa1(2)-p0(2),0, 'k:');
% hold on;
% arm1.plot(q1(1:3)')
% arm2.plot([q1(4:5)',0])
% pause();

% second order forward kinematics
[t0rd2] = fsecondkine5(qdda,qd,q0,Ar,J)
pa2 = p0 +   t0r * dt + t0rd2 * dt^2 / 2;
disp('platform acceleration');
t0rd2
% quiver(p0(1),p0(2),pa2(1)-p0(1),pa2(2)-p0(2),0, 'k:');
% hold on;
% arm1.plot(q1(1:3)')
% arm2.plot([q1(4:5)',0])

% other verifications
Jinv_ = inv(J);
if isequal(Jinv_,Jinv)
    disp('Jinv = inv(J)');
else
        Jinv
        Jinv_
end
l = 0.205;
d1 = [-0.125, l, l];
d2 = [0.125, l, l];
q1 = q0(1:3);
q2 = q0(4:5);
qd1 = qd(1:3);
qd2 = qd(4:5);
qdd1 = qdd(1:3);
qdd2 = qdd(4:5);
% Twists
tw11 = [-d1(2)*sin(q1(1))-d1(3)*sin(q1(1)+q1(2)), d1(2)*cos(q1(1))+d1(3)*cos(q1(1)+q1(2)), 0, 0, 0, 1]';
tw12 = [-d1(3)*sin(q1(1)+q1(2)), d1(3)*cos(q1(1)+q1(2)), 0, 0, 0, 1]';
tw13 = [0, 0, 0, 0, 0, 1]';
tw21 = [-d2(2)*sin(q2(1)) - d2(3)*sin(q2(1)+q2(2)), d2(2)*cos(q2(1)) + d2(3)*cos(q2(1)+q2(2)), 0, 0, 0, 1]';
tw22 = [-d2(3)*sin(q2(1)+q2(2)), d2(3)*cos(q2(1)+q2(2)), 0, 0, 0, 1]';
disp('Verificacion velocidades')
t0p = Psit * t0r
t0p1 = [tw11 tw12 tw13]*qd1
t0p2 = [tw21 tw22]*qd2
disp('Verificacion aceleraciones')
t0pd = Psit * t0rd + Psidt * t0r
t0pd1 = [tw11 tw12 tw13]*qdd1 + b0p1
t0pd2 = [tw21 tw22]*qdd2 + b0p2
