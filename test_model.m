close all;
% q1 = [0.4072    1.5951    0.7321]';
% q2 = [1.1393    1.5951         0]';
% qda = [0.5; 0.5];
% qdda =[1;1];
% n=50;
% % x = -0.1:0.2/n:0.1;
% % y = 0.2:0.1/n:0.3;
% p= zeros(2,1);
% % g = -9.80665;
% g = 0;
% t0r = p;
% t0rd = p;
% % figure(5)
% % lspb(0.2,0.3,n);
% % [p(:,1),t0r(:,1),t0rd(:,1)] = lspb(-0.1,0.1,n);
% % [p(:,2),t0r(:,2),t0rd(:,2)] = lspb(0.35,0.3,n);
% x0 = 0;
% y0 = 0.35;
% p0 = [x0;y0];
% t0r0 = [0;0];
% t0rd0 = [0;0];
%     q_mode = ikine5(p0);
%     q0 = q_mode(:,4);
%     q = q_mode(:,4); % Using -+ mode (4)
%     qa = [q(1);q(4)];
%     % p_ = fkine5(q(:,i))
% %     if q(2,i) == 0 || q(2,i) ==pi || q(2,i) ==-pi || q(2,i) ==-2*pi ...
% %             || q(2,i) ==2*pi || q(5,i)==0 || q(5,i)==pi || q(5,i)==-pi ...
% %             || q(5,i)==-2*pi || q(5,i)==2*pi
% %         disp('Type 1 singularity found');
% %     elseif q(1,i)+q(2,i)-q(4,i)-q(5,i)==0 || q(1,i)+q(2,i)-q(4,i)-q(5,i)==pi ...
% %             || q(1,i)+q(2,i)-q(4,i)-q(5,i)==-pi || q(1,i)+q(2,i)-q(4,i)-q(5,i)==2*pi ...
% %             || q(1,i)+q(2,i)-q(4,i)-q(5,i)==-2*pi
% %         disp('Type 2 singularity found')
% %     end
% %     q(1,i)+q(2,i)-q(4,i)-q(5,i)
%     % t0r = ffirstkine5(q,qda);
%     [qd,qda,~,Ar,B,J,Jinv,Jt,Jta,Jtd,Psit] = ifirstkine5(t0r,p0,q0);
%     % qda_
% 
%     % t0rd = fsecondkine5(qdda,qd,q,Ar,J);
%     [qdd,qdda,qddu,Psit1,Psit2,atp,ad,b0p1,b0p2,Psidt] = isecondkine5(t0rd,t0r,p0,q0,qd,Ar,B,Jinv,Jt,Jta,Jtd,Psit);
    % qdda_
%     T(:,i) = idyn5(q(:,i),qd(:,i),qdd(:,i),t0rd(i,:)',J,Jt,Jta,Jtd,Psit);
%     qdda_d(:,i) = ddyn5(T(:,i),q(:,i),qd(:,i),qdd(:,i),t0rd(i,:)',J,Jt,Jta,Jtd,Psit,atp,ad);
% qa = [q(1,:);q(4,:)];
% x0 = p(1,1);
% y0 = p(1,2);
% planning a circle
dt = 10e-3;
% ang = 2*pi;
% qv = 1;
% qa = 2;
% PSt = [0.15;0.26];
% P0C = [0.15;0.23];
%     
% 
% qaux_ = ikine5(PSt);
% show working modes
% for i=1:4
% arm1.plot(qaux_(1:3,i)')
% hold on;
% arm2.plot([qaux_(4:5,i)',0])
% pause();
% end


% qprev = qaux_(:,4);
% [q, qd, qdd, p, t0r, t0rd, tmax] = cirtraj(PSt ,P0C, qprev, ang, qv, qa, dt);
% v = ones(2,1);
% a = 2 * ones(2,1);
% p1 = [0.1;0.3];
% p2 = [0.0;0.3];
% p3 = [-0.1;0.3];
% 
% q0 = ikine5(p1);
% % q1 = ikine5(p3);
% % q_ = [q0(1,4), q1(1,4); q0(4,4), q1(4,4)];
% qprev = q0(:,4);
% % 
% [q, qd, qdd, p, t0r, t0rd, tmax] = linctraj([p1,p3], qprev, v, a, dt);

b = 10;
x = 0.1:(0.15-0.1)/b:0.15;
y = 0.32:(0.28-0.32)/b:0.28;
for i=1:length(x)
    if rem(i,2) == 0
        q0 = ikine5([x(1);y(i)]);
        q_(:,i) = [q0(1,4);q0(4,4)];
    else
        q0 = ikine5([-x(1);y(i)]);
        q_(:,i) = [q0(1,4);q0(4,4)];
    end
end
% plot(q_(1,:))
qv = ones(length(q_(1,:)),1) * 0.2;
qa = ones(length(q_(1,:)),1) * 0.5;
% q0 = ikine5(p1);
% qprev = q0(:,4);
% q1 = ikine5(p2);
% qprev1 = q1(:,4);
% q2 = ikine5(p3);
% qprev2 = q2(:,4);
% qadesire = [qprev(1),qprev1(1), qprev2(1); qprev(4), qprev1(4), qprev2(4)]; 
[q, qd, qdd, p, t0r, t0rd, tmax] = jotraj(q_, qv, qa, dt);


% [p(:,2),t0r(:,2),t0rd(:,2)] = lspb(0.35,0.3,n);
x0 = p(1,1);
y0 = p(2,1);
p0 = p(:,1);
t0r0 = t0r(:,1);
t0rd0 = t0rd(:,1);
q0 = q(:,1);
qd0 = qd(:,1);
qdd0 = qdd(:,1);
% PID Tunning
aj = 0.0562;
Fv = 1.5e-1;
wpos = Fv / aj;
wj = 15 * wpos;
Kp = 3*aj*wj^2 * eye(2); 
Kd = (3*aj*wj-Fv) * eye(2);
Ki = aj*wj^3 * eye(2);


% PD Inverse dynamic control
l = 0.205;
m11 = 0.126; %kg
mx11 = m11 * 0.1025;
zz11 = 0.001774; %kg-m^2;

m12 = 0.126; %kg
mx12 = m12 * 0.1025;
zz12 = 0.001774; %kg-m^2;

m21 = 0.126; %kg
mx21 = m21 * 0.1025;
zz21 = 0.001774; %kg-m^2;

m22 = 0.126; %kg
mx22 = m22 * 0.1025;
zz22 = 0.001774; %kg-m^2;

Ia11 = 4.8e-5; % kg-m2
Ia21 = 4.8e-5; % kg-m2

m4 = 0.272;

Jeq2 = Ia21 + zz21 + zz12 + l^2*m4 + l^2*m22 + 2 * l * mx22; % cos = 1 and sin = 0
Jeq1 = Ia11 + zz11 + zz22 + l^2*m4 + l^2*m12 + 2 * l * mx12; % ver si es + o -
beq1 = Fv;
beq2 = Fv;


% Kpi = [(2  + 0.75) * Jeq1 * wpos1^2; (2 + 0.75) * Jeq2 * wpos2^2];
% this are working good. 
% wn1 = beq1 / Jeq1;
% wpos1 =   20 * wn1;
% wn2 = beq2 / Jeq2;
% wpos2 =  20 * wn2;
% zita = 1;
% Kpi = [ wpos1^2;  wpos2^2];
% Kdi = [2 * zita *  wpos1; 2  * zita * wpos2];
% Adding integral accion.
wn1 = beq1 / Jeq1;
wpos1 =   15 * wn1;
wn2 = beq2 / Jeq2;
wpos2 =  15 * wn2;
n = 3;
Kii = [wpos1^3; wpos2^3];
Kpi = [n * wpos1^2; n * wpos2^2];
Kdi = [n *  wpos1; n * wpos2];

% figure(1)
% plot(T(1,:),'g--');
% hold on;
% plot(T(2,:),'k-.');
% title('Torque')
% 
% figure(2)
% plot(qdda_d(1,:),'r-.');
% hold on;
% plot(qdda(1,:),'b--');
% title('Acelerations q11')
% 
% figure(3)
% plot(qdda(2,:),'m--');
% hold on;
% plot(qdda_d(2,:),'g-.');
% title('Acelerations q21')
% 
% figure(4)
% plot(q(1,:),'r-.');
% hold on;
% plot(q(4,:),'b--');

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
% for i=1:4
%     if isnan(q_(1,i))
%         continue;
%     end
% Plot circle for verify
% figure(5)
% plot(p(1,:),p(2,:),'r.');
% hold on;
% circle([0,0.3], .03)
% grid on;

for i=1:200:length(simout.Data)
figure(5)
% plot(p(:,1),p(:,2));
% p_(:,i) = fkine5([simout.Data(i,1);simout.Data(i,4)] ,1);
arm1.plot(simout.Data(i,1:3))
hold on;
%trplot(arm1.base * arm1.links(1).A(q(1,1)) )
arm2.plot([simout.Data(i,4:5) 0])
end
% Plot circle for verify 
% plot(p_(1,:),p_(2,:),'g+');
% hold on
% plot(p(1,:),p(2,:),'r.');


% %trplot(arm2.base * arm2.links(1).A(q(4,1)) )

% fkine(arm2,[q(4:5,1)' 0])
% fkine(arm1,q(1:3,1)')
