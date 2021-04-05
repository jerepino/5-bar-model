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
n=50;
x=  zeros(n,1); %-0.1 * ones(n,1);
dt = 0.3-0.1;
y = 0.30:-dt/n:0.10;
j = 1;
for i=1:n
% p0 = [0; 0.19];
p0 = [x(i);y(i)];
q0 = ikine5(p0);
% show working modes
% for i=1:4
% arm1.plot(q0(1:3,i)')
% hold on;
% arm2.plot([q0(4:5,i)',0])
% pause();
% end
q0 = q0(:,4); % Using -+ mode (4)
    l = 0.205;
    l0 = 0.125;
    A12 = [l*cos(q0(1))-l0; l*sin(q0(1))];
    A22 = [l*cos(q0(4))+l0; l*sin(q0(4))];

    if sqrt((A12-A22)' * (A12-A22)) >= (2*l-0.0001)
        q_sing(:,j) = [q0(1);q0(4)];
        index(j) = i;
        j = j+1;
    end
qa(:,i) = [q0(1);q0(4)];

% Forward kinematics validation
% [p1, q1] = fkine5(qa0,1);
% [p2, q2] = fkine5(qa0,-1);

% first order inverse kinematics
t0r = [0;0];
% dt = 0.1;
% pv1 = p0 + t0r * dt;
[qd,qda(:,i),qdu,Ar,B,J,Jinv,Jt,Jta,Jtd,Psit] = ifirstkine5(t0r,p0,q0);

% first order forward kinematics
% [t0r2] = ffirstkine5(q0,qda);
% pv2 = p0 + t0r2 * dt;

% second order inverse kinematics
t0rd = [0;0];
% pa1 = p0 +   t0r * dt + t0rd * dt^2 / 2;
[qdd,~,~,~,~,atp,ad,~,~,~] = isecondkine5(t0rd,t0r,p0,q0,qd,Ar,B,Jinv,Jt,Jta,Jtd,Psit);

% second order forward kinematics
% [t0rd2] = fsecondkine5(qdda,qd,q0,Ar,J)
% pa2 = p0 +   t0r * dt + t0rd2 * dt^2 / 2;

% dynamic
T = [0;0];
[qdda_dynamics(:,i), Ta] = ddyn5(T,q0,qd,qdd,t0rd,J,Jt,Jta,Jtd,Psit,atp,ad);
% Ta
% [t0rd2] = fsecondkine5(qdda_dynamics(:,i),qd,q0,Ar,J);
t0rd2 = [0;-9.80665];
[qdd,qdda,~,~,~,~,~,~,~,~] = isecondkine5(t0rd2,t0r,p0,q0,qd,Ar,B,Jinv,Jt,Jta,Jtd,Psit);
[T_(:,i)] = idyn5(q0,qd,qdd,t0rd,J,Jt,Jta,Jtd,Psit);

g = -9.80665;
t0rd1 = [0;g];
[qdd,qdda_kinematics(:,i),qddu,Psit1,Psit2,atp,ad,b0p1,b0p2,Psidt] = isecondkine5(t0rd1,t0r,p0,q0,qd,Ar,B,Jinv,Jt,Jta,Jtd,Psit);
% qdda_kinematics
end

T = zeros(2,length(T_));
%Integrations 
qda_dynamics(1,:) = cumtrapz(qdda_dynamics(1,:));
qda_kinematics(1,:) = cumtrapz(qdda_kinematics(1,:));
qa_dynamics(1,:) = cumtrapz(qda_dynamics(1,:));
qa_kinematics(1,:) = cumtrapz(qda_kinematics(1,:));

qda_dynamics(2,:) = cumtrapz(qdda_dynamics(2,:));
qda_kinematics(2,:) = cumtrapz(qdda_kinematics(2,:));
qa_dynamics(2,:) = cumtrapz(qda_dynamics(2,:));
qa_kinematics(2,:) = cumtrapz(qda_kinematics(2,:));

% Plots
f1 = figure(1);
p1 = uipanel('Parent',f1,'BorderType','none'); 
p1.Title = 'Accelerations'; 
p1.TitlePosition = 'centertop'; 
subplot(2,1,1,'Parent',p1)
plot(qdda_dynamics(1,:));
hold on;
grid on;
plot(qdda_kinematics(1,:));
plot(index,q_sing(1,:),'mo');
plot(index,q_sing(2,:),'mo');
hl1 = legend('$\ddot{qa}_{11dynamics}$','$\ddot{qa}_{11kinematics}$');
set(hl1, 'Interpreter', 'latex');
title('Qdda11');
subplot(2,1,2,'Parent',p1)
plot(qdda_dynamics(2,:));
hold on;
grid on;
plot(qdda_kinematics(2,:));
plot(index,q_sing(1,:),'mo');
plot(index,q_sing(2,:),'mo');
hl2 = legend('$\ddot{qa}_{21dynamics}$','$\ddot{qa}_{21kinematics}$');
set(hl2, 'Interpreter', 'latex');
set(hl2, 'Location','southeast');
title('Qdda21');
 

f2 = figure(2);
p2 = uipanel('Parent',f2,'BorderType','none'); 
p2.Title = 'Velocities'; 
p2.TitlePosition = 'centertop'; 
subplot(3,1,1,'Parent',p2)
plot(qda_dynamics(1,:));
hold on;
grid on;
plot(qda_kinematics(1,:));
hl1 = legend('$\dot{qa}_{11dynamics}$','$\dot{qa}_{11kinematics}$');
set(hl1, 'Interpreter', 'latex');
set(hl1, 'Location','southeast');
title('Qda11');
subplot(3,1,2,'Parent',p2)
plot(qda_dynamics(2,:));
hold on;
grid on;
plot(qda_kinematics(2,:));
hl2 = legend('$\dot{qa}_{21dynamics}$','$\dot{qa}_{21kinematics}$');
set(hl2, 'Interpreter', 'latex');
title('Qda21');
subplot(3,1,3,'Parent',p2)
plot(qda(1,:));
hold on;
grid on;
plot(qda(2,:));
hl2 = legend('$\dot{qa}_{11}$','$\dot{qa}_{21}$');
set(hl2, 'Interpreter', 'latex');
set(hl2, 'Location','southeast');
title('ikinematics');
 

f3 = figure(3);
p3 = uipanel('Parent',f3,'BorderType','none'); 
p3.Title = 'Positions'; 
p3.TitlePosition = 'centertop'; 
subplot(3,1,1,'Parent',p3)
plot(qa_dynamics(1,:));
hold on;
grid on;
plot(qa_kinematics(1,:));
hl1 = legend('$qa_{11dynamics}$','$qa_{11kinematics}$');
set(hl1, 'Interpreter', 'latex');
set(hl1, 'Location','southeast');
title('Qa11');
subplot(3,1,2,'Parent',p3)
plot(qa_dynamics(2,:));
hold on;
grid on;
plot(qa_kinematics(2,:));
hl2 = legend('$qa_{21dynamics}$','$qa_{21kinematics}$');
set(hl2, 'Interpreter', 'latex');
title('Qa21');
subplot(3,1,3,'Parent',p3)
plot(qa(1,:));
hold on;
grid on;
plot(qa(2,:));
hl2 = legend('$qa_{11}$','$qa_{21}$');
set(hl2, 'Interpreter', 'latex');
plot(index,q_sing(1,:),'mo');
plot(index,q_sing(2,:),'mo');
% set(hl2, 'Location','southeast');
title('ikinematic');

f4 = figure(4);
p4 = uipanel('Parent',f4,'BorderType','none'); 
p4.Title = 'Torque'; 
p4.TitlePosition = 'centertop'; 
subplot(2,1,1,'Parent',p4)
plot(T(1,:));
hold on;
grid on;
plot(T_(1,:));
hl1 = legend('T1','$T1_{1}$');
plot(index,q_sing(1,:),'mo');
plot(index,q_sing(2,:),'mo');
set(hl1, 'Interpreter', 'latex');
set(hl1, 'Location','southeast');
title('T1');
subplot(2,1,2,'Parent',p4)
plot(T(2,:));
hold on;
grid on;
plot(T_(2,:));
hl2 = legend('T2','$T2_{2}$');
plot(index,q_sing(1,:),'mo');
plot(index,q_sing(2,:),'mo');
set(hl2, 'Interpreter', 'latex');
title('T2');
