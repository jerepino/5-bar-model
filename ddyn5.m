function [qdda] = ddyn5(T,q,qd,qdd,t0rd,J,Jt,Jta,Jtd,Psit,atp,ad)

    l = 0.205;
    d1 = [-0.125, l, l];
    d2 = [0.125, l, l];
    q1 = q(1:3);
    q2 = q(4:5);
    qd1 = qd(1:3);
    qd2 = qd(4:5);
    qdd1 = qdd(1:3);
    qdd2 = qdd(4:5);
    xdd = t0rd(1);
    ydd = t0rd(2);
    % Numerical Solution Dynamics
    % Gravity 
    g = 9.80665;
    
    m11 = 83.731e-3; %kg
    mx11 = m11 * 8.82e-2;
    my11 =  0; %m11 * 1.44e-2;
    zz11 = 12375.675e-7;
    
    m21 = 102.422e-3; %kg
    mx21= m21 * 7.21e-2;
    my21= 0; %m21 * 1.48e-2;
    zz21 = 12469.254e-7;
    
    m12 = 72.02e-3;
    my12 = 0;
    mx12 = m12 * 10.25e-2;
    zz12 = 12283.918e-7;
    
    m22 = 72.02e-3;
    my22 = 0;
    mx22 = m22 * 10.25e-2;
    zz22 = 12283.918e-7;
    
    %Fricciones fv = viscosa fs= seca
    fv1 = [1.5e-5 0 0];
    fs1 = [1.5e-5 0 0];
    fv2 = [0 0];
    fs2 = [0 0];
    
%     zz11 = 12375.675e-7; % kg-m²
    Ia11 = 0.48 * 0.0001; % kg-m2
    zz11r = zz11 + Ia11 + d1(2)^2 * m12;
%     zz11r = 2.11e-2;

%     zz21 = 12469.254e-7; % kg-m²
    Ia21 = 0.48 * 0.0001; % kg-m2
    zz21r = zz21 + Ia21 + d2(2)^2 * m22;
%     zz21r = 2.24e-2;
    m4 = 0.272;
    xdd = t0rd(1);
    ydd = t0rd(2);
    % Direct dynamics model
    M11 = zz11r + zz12 + 2 * d1(2) * mx12 * cos(q1(2)) + ...
          2 * d1(2) * my12 * sin(q1(2));

    M12 = zz12 + d1(2) * mx12 * cos(q1(2)) + ...
          d1(2) * my12 * sin(q1(2));
    M21 = M12;

    M44 = zz21r + zz22 + 2 * d2(2) * mx22 * cos(q2(2)) + ...
          2 * d2(2) * my22 * sin(q2(2));

    M45 = zz22 + d2(2) * mx22 * cos(q2(2)) + ...
          d2(2) * my22 * sin(q2(2));
    M54 = M45;

    Mt = [M11 M12  0 0     0;
          M21 zz12 0 0     0;
          0   0    0 0     0;
          0   0    0 M44 M45;
          0   0    0 M54 zz22];

    c1 = -d1(2) * mx12 * qd1(2) * (2 * qd1(1) + qd1(2)) * sin(q1(2)) + ...
         d1(2) * my12 * qd1(2) * (2 * qd1(1) + qd1(2)) * cos(q1(2)) + ...
         fs1(1) * sign(qd1(1)) + fv1(1) * qd1(1) + ...
         (my12*sin(q1(1) + q1(2)) - mx12*cos(q1(1) + q1(2)) - mx11*cos(q1(1)) + ...
         my11*sin(q1(1)) - d1(2)*m12*cos(q1(1)))*g;

    c2 = d1(2) * mx12 * qd1(1)^2 * sin(q1(2)) - ...
         d1(2) * my12 * qd1(1)^2 * cos(q1(2)) + ...
         fs1(2) * sign(qd1(2)) + fv1(2) * qd1(2) + ...
         (my12*sin(q1(1) + q1(2)) - mx12*cos(q1(1) + q1(2)))*g;

    c3 = fs1(3)*sign(qd1(3)) + fv1(3)*qd1(3);

    c4 = -d2(2) * mx22 * qd2(2) * (2 * qd2(1) + qd2(2) ) * sin(q2(2)) + ...
         d2(2) * my22 * qd2(2) * (2 * qd2(1) + qd2(2) ) * cos(q2(2)) + ...
         fs2(1) * sign( qd2(1) ) + fv2(1) * qd2(1) + ...
         (my22*sin(q2(1) + q2(2)) - mx22*cos(q2(1) + q2(2)) - mx21*cos(q2(1)) ...
         + my21*sin(q2(1)) - d2(2)*m22*cos(q2(1)))*g;

    c5 = d2(2) * mx22 * qd2(1)^2 * sin(q2(2)) - ...
         d2(2) * my22 * qd2(1)^2 * cos(q2(2)) + ...
         fs2(2) * sign(qd2(2)) + fv2(2) * qd2(2) + ...
         (my22*sin(q2(1) + q2(2)) - mx22*cos(q2(1) + q2(2)))*g;

    ct = [c1; c2; c3; c4; c5];

    Et = [1 0 0 0 0;
          0 0 1 0 0;
          0 0 0 1 0;
          0 1 0 0 0;
          0 0 0 0 1];

    M0p = m4 * [eye(2), zeros(2,4); zeros(4,2), zeros(4)];
    % c0p = zeros(6,1);
    w0p_ = M0p * [xdd; ydd-g; zeros(4,1)]; % + c0p;

    % M0r = Psit' * M0p * Psit; % Igual a m4*eye(2)
    M0r = m4 * eye(2);

    Jd = inv(Jtd) * (Jt * J - Jta);

    M = [eye(2), Jd'] * Et'*Mt*Et * [eye(2);Jd] + J'*M0r*J;

        % wr1 = [cos(q1(1)+q1(2)); sin(q1(1)+q1(2)); 0; 0; 0; 0];
        % wr2 = [cos(q2(1)+q2(2)); sin(q2(1)+q2(2)); 0; 0; 0; 0];

        % b0p1 = [-d1(2)*qd1(1)*cos(q1(1)) - d1(3)*cos(q1(1)+q1(2))*(qd1(1)+qd1(2)); -d1(2)*qd1(1)*sin(q1(1)) - d1(3)*sin(q1(1)+q1(2))*(qd1(1)+qd1(2)); 0; 0; 0; 0] * qd1(1) + [- d1(3)*cos(q1(1)+q1(2));- d1(3)*sin(q1(1)+q1(2)); 0; 0; 0; 0] * qd1(2)*(qd1(1)+qd1(2));
        % b0p2 = [-d2(2)*qd2(1)*cos(q2(1)) - d2(3)*cos(q2(1)+q2(2))*(qd2(1)+qd2(2)); -d2(2)*qd2(1)*sin(q2(1)) - d2(3)*sin(q2(1)+q2(2))*(qd2(1)+qd2(2)); 0; 0; 0; 0] * qd2(1) + [- d2(3)*cos(q2(1)+q2(2));- d2(3)*sin(q2(1)+q2(2)); 0; 0; 0; 0] * qd2(2)*(qd2(1)+qd2(2));

        % b0p = [wr1'*b0p1; wr2'*b0p2];
    % b0pr = b0p; % - A * Psidt * t0rn; % Its b0p

    % at_ = -inv(Ar)*b0pr;
    % atp = Psit * at + Psidt * t0rn;
    % atp_ = Psit * at_ + Psidt * t0rn

    c0r = M0p * atp; % Debe ser igual a 0
    % c0r = (M0p * atp_) % Debe ser igual a 0
    % c0r = zeros(6,1);

    % ad = inv(Jtd) * (dc + Jt*at);
    % ad_ = inv(Jtd) * (dc + Jt * at_)

    c = [eye(2) Jd'] * Et' * (ct + Mt*Et*[zeros(2,1); ad]) + J'*Psit'*c0r; % if c0r and c0p are 0 the second term is 0

    % TT = M * qdda + c
    % qdda_ = inv(M) * (TT-c)
    % T
    qdda = inv(M) * (T-c);




