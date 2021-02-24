function [T] = idyn5(q,qd,qdd,t0rd,J,Jt,Jta,Jtd,Psit)

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
    % Dynamic
    % Numerical Solution Dynamic
    m12 = 83.731e-3; %kg
    mx12 = m12 * 8.82e-2;
    my12 = m12 * 1.44e-2;
    m22 = 102.422e-3; %kg
    mx22= m22 * 7.21e-2;
    my22= m22 * 1.48e-2;

    zz12 = 12283.918e-7;
    zz22 = 12283.918e-7;
    %Fricciones fv = viscosa fs= seca
    fv1 = [1.5e-5 0 0];
    fs1 = [1.5e-5 0 0];
    fv2 = [0 0];
    fs2 = [0 0];
    
    zz11 = 12375.675e-7; % kg-m²
    Ia11 = 0.48 * 0.0001; % kg-m2
    zz11r = zz11 + Ia11 + d1(2)^2 * m12;
%     zz11r = 2.11e-2;

    zz21 = 12469.254e-7; % kg-m²
    Ia21 = 0.48 * 0.0001; % kg-m2
    zz21r = zz21 + Ia21 + d2(2)^2 * m22;
%     zz21r = 2.24e-2;
    m4 = 0.272;

    Tt1(1) = zz11r * qdd1(1) + zz12 * (qdd1(1) + qdd1(2)) + ...
        d1(2) * mx12 * ((2 * qdd1(1) + qdd1(2)) * cos(q1(2)) - qd1(2) * (2 * qd1(1) + qd1(2)) * sin(q1(2))) + ...
        d1(2) * my12 * ((2 * qdd1(1) + qdd1(2)) * sin(q1(2)) + qd1(2) * (2 * qd1(1) + qd1(2)) * cos(q1(2))) + ...
        fs1(1) * sign(qd1(1)) + fv1(1) * qd1(1);
    Tt1(2) = zz12 *(qdd1(1) + qdd1(2)) + d1(2)*mx12*(qdd1(1)*cos(q1(2))+ qd1(1)^2 * sin(q1(2))) + ...
        d1(2)*my12*(qdd1(1)*sin(q1(2))- qd1(1)^2*cos(q1(2))) + fs1(2)*sign(qd1(2)) + fv1(2)*qd1(2);
    Tt1(3) = fs1(3) * sign(qd1(3)) + fv1(3)*qd1(3);

    Tt2(1) = zz21r * qdd2(1) + zz22 * (qdd2(1) + qdd2(2)) + ...
             d2(2) * mx22 * ((2 * qdd2(1) + qdd2(2)) * cos(q2(2)) - qd2(2) * (2 * qd2(1) + qd2(2)) * sin(q2(2))) + ...
             d2(2) * my22 * ((2 * qdd2(1) + qdd2(2)) * sin(q2(2)) + qd2(2) * (2 * qd2(1) + qd2(2)) * cos(q2(2))) + ...
             fs2(1) * sign(qd2(1)) + fv2(1) * qd2(1);
    Tt2(2) = zz22 *(qdd2(1) + qdd2(2)) + d2(2) * mx22 * (qdd2(1) * cos(q2(2)) + qd2(1)^2 * sin(q2(2))) + ...
             d2(2) * my22 * (qdd2(1) * sin(q2(2)) - qd2(1)^2 * cos(q2(2))) + fs2(2) * sign(qd2(2)) + fv2(2) * qd2(2);
    % End effector like a punctual mass fixed on B22
    w0p = m4 * [xdd; ydd; zeros(4,1)];
    w0r = Psit' * w0p;
    % IDM
    Tta = [Tt1(1); Tt2(1)];
    Ttd = [Tt1(2:3)'; Tt2(2)];

    Jd = inv(Jtd) * (Jt * J - Jta);

    T = Tta + J' * w0r + Jd' * Ttd;
