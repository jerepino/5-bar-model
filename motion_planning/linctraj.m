function [q, qd, qdd, p, t0r, t0rd, tmax] = linctraj(pd, qprev, v, a, dt)
%% Funcion para generar trayectorias lineales entre puntos cartesianos
% pd = matriz 2XI con las posiciones deseadas
% v = maxima velocidad cartesiana
% a = maxima aceleracion cartesiana
% dt = delta time
%%
n = length(pd(1,:));

    
    if nargin < 3
        for i=1:n
            v(i) = 1;
            a(i) = 0.5;
        end
        dt = 0.05;
    elseif nargin < 5
        dt = 0.05;
    end

    T = zeros(1, n-1);
    tau = zeros(1, n-1);

    for j = 1:n-1
        dis_p = sqrt(sum((pd(:,j+1) - pd(:,j)).^2));
    %     dis_fi = sqrt(sum((fid(j+1,:) - fid(j,:)).^2));
        dq = max(dis_p); %dq = desplazamiento total
        [T(j),tau(j)] = tlparam(dq, v(j), a(j), dt);
    end

    tmax = sum(T) + tau(end)+dt;
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
    for t = 0 : dt : tmax
        s = sfun(t - T_ant, T(k), tau(k));
        ds = sdfun(t - T_ant, T(k), tau(k));
        dds = sddfun(t - T_ant, T(k), tau(k));
        if n > 2 
            T2 = sum(T(1:k));
            s2 = sfun(t-T2, T(k+1), tau(k+1));
            ds2 = sdfun(t-T2, T(k+1), tau(k+1));
            dds2 = sddfun(t-T2, T(k+1), tau(k+1));
            pdiff2 = pd(:,k+2) - pd(:,k+1);
        else
            s2 = 0;
            ds2 = 0;
            dds2 = 0;
            pdiff2 = 0;
        end
        pdiff1 = pd(:,k+1) - pd(:,k);
        p(:,l) = pd(:,k) + s * pdiff1 + s2 * pdiff2;
        t0r(:,l) = ds * pdiff1 + ds2 * pdiff2;
        t0rd(:,l) = dds * pdiff1 + dds2 * pdiff2;

        [q(:,l),~,~] = ikine5(p(:,l),qprev);
        [qd(:,l),~,~,Ar,B,~,Jinv,Jt,Jta,Jtd,Psit] = ifirstkine5(t0r(:,l),p(:,l), ...
                                                                q(:,l));
        [qdd(:,l),~,~,~,~,~,~,~,~,~] = isecondkine5(t0rd(:,l),t0r(:,l),p(:,l), ...
                                        q(:,l),qd(:,l),Ar,B,Jinv,Jt,Jta,Jtd,Psit);
        qprev = q(:,l);

        
        if n > 2 && (s==1) && (k ~= n-2)
            T_ant = T2;
            k = k+1;
        end
        l = l + 1;
    end