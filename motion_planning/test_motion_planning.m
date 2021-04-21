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
qv = ones(3,1);
qa = ones(3,1) * 2;
PSt = [0;0.35];
P0C = [0;0.3];
    
p1 = [-0.0;0.39];
p2 = [-0.01;0.38];
p3 = [-0.1;0.38];

b = 50;
x = 0:(0.15-0)/50:0.15;
y = 0.39:(0.25-0.39)/b:0.25;
for i=1:length(x)
    if rem(i,2) == 0
        q0 = ikine5([x(i);y(i)]);
        q_(:,i) = [q0(1,4);q0(4,4)];
    else
        q0 = ikine5([-x(i);y(i)]);
        q_(:,i) = [q0(1,4);q0(4,4)];
    end
end
qv = ones(length(q_),1);
qa = ones(length(q_),1) * 2;
% q0 = ikine5(p1);
% qprev = q0(:,4);
% q1 = ikine5(p2);
% qprev1 = q1(:,4);
% q2 = ikine5(p3);
% qprev2 = q2(:,4);
% qadesire = [qprev(1),qprev1(1), qprev2(1); qprev(4), qprev1(4), qprev2(4)]; 
[q, qd, qdd, p, t0r, t0rd, tmax] = jotraj(q_, qv, qa, dt);
% [q, qd, qdd, p, t0r, t0rd, tmax] = linctraj([p1,p2,p3], qprev, qv, qa, dt);

% [q, qd, qdd, p, t0r, t0rd, tmax] = cirtraj(PSt ,P0C, qprev, ang, 1, 2, dt);
%% Find max inertia
% for i = 1:length(q(1,:))
%     [~,qda,~,Ar,B,J,Jinv,Jt,Jta,Jtd,Psit] = ifirstkine5(t0r(:,i),p(:,i),q(:,i));
%         % qda_
% 
%         % t0rd = fsecondkine5(qdda,qd,q,Ar,J);
%     [~,qdda,qddu,Psit1,Psit2,atp,ad,b0p1,b0p2,Psidt] = isecondkine5(t0rd(:,i),t0r(:,i),p(:,i),q(:,i),qd(:,i),Ar,B,Jinv,Jt,Jta,Jtd,Psit);
%     [qdda_d(:,i), ~,M{i}] = ddyn5([0;0],q(:,i),qd(:,i),qdd(:,i),t0rd(:,i),J,Jt,Jta,Jtd,Psit,atp,ad);
% end
% 
% for i = 1:length(q(1,:))
%     MM = M{i}
%     M1(i) = max(MM(1,:));
%     M2(i) = max(MM(2,:));
% end
% max(M1)
% max(M2)

%% Plot
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
for i=1:5:length(q)
arm1.plot(q(1:3,i)');
hold on;
arm2.plot([q(4:5,i)' 0]);
end