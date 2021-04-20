function [q, qd, qdd, p, t0r, t0rd, tmax] = jotraj(qadesire, qvel, qacc, dt)
%% Funcion para generar trayectorias entre puntos articulares
% qd = vector de posiciones articulares
% qv = maxima velocidad articular por trayecto
% qa = maxima aceleracion articular por trayecto
%%
n = length(qadesire(1,:));

if nargin < 2
    % velicidades articulares maximas
    dqlim = [deg2rad(360) deg2rad(360) deg2rad(360) deg2rad(360) deg2rad(360)];
    ddqlim = dqlim / 0.8; % Se asume que el robot llegara en 0.5 s a su velocidad maxima linealmente -> a = dqlim/0.5
    for i=1:n
        qvel(i) = min(dqlim);
        qacc(i) = min(ddqlim);
    end
    dt = 0.05;
elseif nargin < 4
    dt = 0.05;
end


T = zeros(1,n-1);
tau = T;

for j = 1:n-1
    dq = max(abs(qadesire(:,j+1) - qadesire(:,j))); %dq = desplazamiento total    
    [T(j),tau(j)] = tlparam(dq, qvel(j), qacc(j), dt);
end
    tmax = sum(T) + tau(end);
    T_ant = 0;
    k = 1;
    l = 1;
    m = round(tmax/dt) + 1;
    p = zeros(2,m);
    t0r = p;
    t0rd = p;
    % T_ = zeros(4,4,m);
    q = zeros(5,m);
    qd = q;
    qdd = q;
    qa = zeros(2,m);
    qda = qa;
    qdda = qa;
    for t = 0 : dt : tmax

        s = sfun(t - T_ant, T(k), tau(k));
        ds = sdfun(t - T_ant, T(k), tau(k));
        dds = sddfun(t - T_ant, T(k), tau(k));
%         dds_(l) = dds;
        if n > 2 
            T2 = sum(T(1:k));
            s2 = sfun(t-T2, T(k+1), tau(k+1));
            ds2 = sdfun(t-T2, T(k+1), tau(k+1));
            dds2 = sddfun(t-T2, T(k+1), tau(k+1));
            qadiff2 = qadesire(:,k+2) - qadesire(:,k+1);
        else
            s2 = 0;
            ds2 = 0;
            dds2 = 0;
            qadiff2 = 0;
        end
        qadiff1 = qadesire(:,k+1) - qadesire(:,k);
        qa(:,l) = qadesire(:,k) + s * qadiff1 + s2 * qadiff2;
        qda(:,l) = ds * qadiff1 + ds2 * qadiff2;
        qdda(:,l) = dds * qadiff1 + dds2 * qadiff2;
        
        [p(:,l),q(:,l)] = fkine5(qa(:,l), 1);
        [t0r(:,l)] = ffirstkine5(q(:,l),qda(:,l));
        [qd(:,l),~,~,Ar,B,J,Jinv,Jt,Jta,Jtd,Psit] = ifirstkine5(t0r(:,l),p(:,l), ...
                                                        q(:,l));
        [t0rd(:,l)] = fsecondkine5(qdda(:,l),qd(:,l),q(:,l),Ar,J);
%         [q(:,l),~,~] = ikine5(p(:,l),qprev);
        [qdd(:,l),~,~,~,~,~,~,~,~,~] = isecondkine5(t0rd(:,l),t0r(:,l),p(:,l), ...
                                        q(:,l),qd(:,l),Ar,B,Jinv,Jt,Jta,Jtd,Psit);
        
%         qprev = q(:,l);

        
        if n > 2 && (s==1) && (k ~= n-2)
            T_ant = T2;
            k = k+1;
        end
        l = l + 1;
    end

%     figure(4)
%     plot(qda(1,:))
%     hold on;
%     plot(qda(2,:))
%     pause()