function [q, qd, qdd, p, t0r, t0rd, tmax] = cirtraj(PSt ,P0C, qprev, ang, qv, qa, dt)
%% Funcion para generar circulos en el plano
% P0C = [xc;yc] distancia desde la base al centro del circulo
% PSt = [xst; yst] posición inicial del circulo
%
%%

if nargin < 4
    qv = 1;
    qa = 0.5;
    dt = 0.01;
    ang = 2 * pi;
elseif nargin < 5
    qv = 1;
    qa = 0.5;
    dt = 0.01;
elseif nargin < 7
    dt = 0.01;    
end

PCst = PSt - P0C;
ro = norm(PCst);
% ro = sqrt(PCst' * PCst)
dq = ang * ro;

[T,tau] = tlparam(dq, qv, qa , dt);
tmax = T+tau+dt;
t = 0 : dt : tmax;

n = length(t);

p = zeros(2,n);
t0r = p;
t0rd = p;

q = zeros(5, n);
qd = zeros(5,n);
qdd = zeros(5,n);

theta0 = atan2(PCst(2), PCst(1));
theta1 = ang + theta0;

for i=1:n
    s = sfun(t(i), T, tau);
    ds = sdfun(t(i), T, tau);
    dds = sddfun(t(i), T, tau);
    
    arcs = theta0 * ro + s * (theta1 - theta0) * ro;
    
    ps = [ro * cos(arcs/ro); ro * sin(arcs/ro)];
    dps = [ -sin(arcs/ro); cos(arcs/ro)];
    ddps = [-cos(arcs/ro)/ro; -sin(arcs/ro)/ro];
    
    p(:,i) = P0C + ps;
    t0r(:,i) = ds * dps;
    t0rd(:,i) = dds * dps + ds^2 * ddps;
    
    [q(:,i),~,~] = ikine5(p(:,i),qprev);
    [qd(:,i),~,~,Ar,B,~,Jinv,Jt,Jta,Jtd,Psit] = ifirstkine5(t0r(:,i),p(:,i), ...
                                                            q(:,i));
    [qdd(:,i),~,~,~,~,~,~,~,~,~] = isecondkine5(t0rd(:,i),t0r(:,i),p(:,i), ...
                                    q(:,i),qd(:,i),Ar,B,Jinv,Jt,Jta,Jtd,Psit);
    qprev = q(:,i);
end