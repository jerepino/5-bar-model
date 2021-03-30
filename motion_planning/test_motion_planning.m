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

dt = 0.01;
ang = 2*pi;
qv = 1;
qa = 0.5;
PSt = [0;0.35];
P0C = [0;0.3];
    

q0 = ikine5(PSt);
qprev = q0(:,4);

[q, qd, qdd, p, t0r, t0rd, tmax] = cirtraj(PSt ,P0C, qprev, ang, qv, qa, dt);

f1 = figure(1);
p1 = uipanel('Parent',f1,'BorderType','none'); 
p1.Title = 'Joint Space'; 
p1.TitlePosition = 'centertop'; 
subplot(3,1,1,'Parent',p1)
plot(qdd(1,:));
hold on;
grid on;
plot(qdd(4,:));
hl1 = legend('$\ddot{qa}_{11}$','$\ddot{qa}_{21}$');
set(hl1, 'Interpreter', 'latex');
title('Acceleration');
subplot(3,1,2,'Parent',p1)
plot(qd(1,:));
hold on;
grid on;
plot(qd(4,:));
hl1 = legend('$\dot{qa}_{11}$','$\dot{qa}_{21}$');
set(hl1, 'Interpreter', 'latex');
title('Velocity');
subplot(3,1,3,'Parent',p1)
plot(q(1,:));
hold on;
grid on;
plot(q(4,:));
hl1 = legend('$qa_{11}$','$qa_{21}$');
set(hl1, 'Interpreter', 'latex');
title('Position');

f2 = figure(2);
p2 = uipanel('Parent',f2,'BorderType','none'); 
p2.Title = 'Task Space'; 
p2.TitlePosition = 'centertop'; 
subplot(3,1,1,'Parent',p2)
plot(t0rd(1,:));
hold on;
grid on;
plot(t0rd(2,:));
hl1 = legend('$\ddot{X}$','$\ddot{Y}$');
set(hl1, 'Interpreter', 'latex');
title('Acceleration');
subplot(3,1,2,'Parent',p2)
plot(t0r(1,:));
hold on;
grid on;
plot(t0r(2,:));
hl1 = legend('$\dot{X}$','$\dot{Y}$');
set(hl1, 'Interpreter', 'latex');
title('Velocity');
subplot(3,1,3,'Parent',p2)
plot(p(1,:));
hold on;
grid on;
plot(p(2,:));
hl1 = legend('$X$','$Y$');
set(hl1, 'Interpreter', 'latex');
title('Position');

figure(3)
plot(p(1,:),p(2,:),'r.');
for i=1:10:length(q)
arm1.plot(q(1:3,i)');
hold on;
arm2.plot([q(4:5,i)' 0]);
end