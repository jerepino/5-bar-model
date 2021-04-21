show_trajectory = false;
dt = 10e-3;
dt_controller = 1e-3;

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
qv = ones(length(q_(1,:)),1) * 100;
qa = ones(length(q_(1,:)),1) * 100;
 
[q, qd, qdd, p, t0r, t0rd, tmax] = jotraj(q_, qv, qa, dt);

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

% this are working good. 
wn1 = beq1 / Jeq1;
wpos1 =   25 * wn1;
wn2 = beq2 / Jeq2;
wpos2 =  25 * wn2;
zita = 1;
Kpi = [ wpos1^2;  wpos2^2];
Kdi = [2 * zita *  wpos1; 2  * zita * wpos2];
% Adding integral accion.
% wn1 = beq1 / Jeq1;
% wpos1 =   15 * wn1;
% wn2 = beq2 / Jeq2;
% wpos2 =  15 * wn2;
% n = 3;
% Kii = [wpos1^3; wpos2^3];
% Kpi = [n * wpos1^2; n * wpos2^2];
% Kdi = [n *  wpos1; n * wpos2];

% Leaky integrator filter 
lambda = 0.7;
% fs = 1/1e-3;
% fs = 30 ;%* wn1 / (2*pi)
step = 2 * pi / 4096; 
Nr = 4096;
% fs  = max(qd(1,:))  / (2*pi) * Nr ; % 314 rad / s tipica velocidad maxima de un motor
fs  = (1/0.1);
fc0 = - log(lambda) / pi * fs;
wco = 2 * pi * fc0;

if show_trajectory
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

    for i=1:200:length(simout.Data)
    figure(5)
    % plot(p(:,1),p(:,2));
    % p_(:,i) = fkine5([simout.Data(i,1);simout.Data(i,4)] ,1);
    arm1.plot(simout.Data(i,1:3))
    hold on;
    %trplot(arm1.base * arm1.links(1).A(q(1,1)) )
    arm2.plot([simout.Data(i,4:5) 0])
    end
end
