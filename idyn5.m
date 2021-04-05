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
    % Gravity 
    g = -9.80665;
    
   % Dynamics parameters
    m11 = 0.126; %kg
    mx11 = m11 * 0.1025;
    my11 =  0; %m11 * 1.44e-2;
    zz11 = 0.001774; %kg-m^2;
    
    m12 = 0.126; %kg
    mx12 = m12 * 0.1025;
    my12 =  0; %m11 * 1.44e-2;
    zz12 = 0.001774; %kg-m^2;
    
    m21 = 0.126; %kg
    mx21 = m21 * 0.1025;
    my21 =  0; %m11 * 1.44e-2;
    zz21 = 0.001774; %kg-m^2;
    
    m22 = 0.126; %kg
    mx22 = m22 * 0.1025;
    my22 =  0; %m11 * 1.44e-2;
    zz22 = 0.001774; %kg-m^2;
    
    Ia11 = 4.8e-5; % kg-m2
%     zz11r = zz11 + Ia11 + d1(2)^2 * m12;
    
    Ia21 = 4.8e-5; % kg-m2
%     zz21r = zz21 + Ia21 + d2(2)^2 * m22;
   
    %Fricciones fv = viscosa fs= seca
    fv1 = [1.5e-1 0 0];
    fs1 = [1.5e-1 0 0];
    fv2 = [0 0];
    fs2 = [0 0];
    % End effector mass
    m4 = 0.272;

    % Inverse dynamics model
    Tt1(1) = (my12*sin(q1(1) + q1(2)) - mx12*cos(q1(1) + q1(2)) ...
        - mx11*cos(q1(1)) + my11*sin(q1(1)) - d1(2)*m12*cos(q1(1)))*g ...
        + Ia11*qdd1(1) + fv1(1)*qd1(1) + fs1(1)*sign(qd1(1)) ...
        + mx12*(- d1(2)*sin(q1(2))*qd1(2)^2 - 2*d1(2)*qd1(1)*sin(q1(2))*qd1(2) ...
        + 2*d1(2)*qdd1(1)*cos(q1(2)) + d1(2)*qdd1(2)*cos(q1(2))) ...
        - my12*(d1(2)*cos(q1(2))*qd1(2)^2 + 2*d1(2)*qd1(1)*cos(q1(2))*qd1(2) ...
        + 2*d1(2)*qdd1(1)*sin(q1(2)) + d1(2)*qdd1(2)*sin(q1(2))) ...
        + qdd1(1)*zz11 + qdd1(1)*zz12 + qdd1(2)*zz12 + m12*qdd1(1)*d1(2)^2;
   
    Tt1(2) = (my12*sin(q1(1) + q1(2)) - mx12*cos(q1(1) + q1(2)))*g ...
        + zz12*(qdd1(1) + qdd1(2)) + fv1(2)*qd1(2) + fs1(2)*sign(qd1(2)) ...
        + mx12*(d1(2)*sin(q1(2))*qd1(1)^2 + d1(2)*qdd1(1)*cos(q1(2))) ...
        - my12*(-d1(2)*cos(q1(2))*qd1(1)^2 + d1(2)*qdd1(1)*sin(q1(2)));

    Tt1(3) = fs1(3) * sign(qd1(3)) + fv1(3)*qd1(3);
   
    Tt2(1) = (my22*sin(q2(1) + q2(2)) - mx22*cos(q2(1) + q2(2)) ...
        - mx21*cos(q2(1)) + my21*sin(q2(1)) - d2(2)*m22*cos(q2(1)))*g ...
        + Ia21*qdd2(1) + fv2(1)*qd2(1) + fs2(1)*sign(qd2(1)) ...
        + mx22*(- d2(2)*sin(q2(2))*qd2(2)^2 - 2*d2(2)*qd2(1)*sin(q2(2))*qd2(2) ...
        + 2*d2(2)*qdd2(1)*cos(q2(2)) + d2(2)*qdd2(2)*cos(q2(2))) ...
        - my22*(d2(2)*cos(q2(2))*qd2(2)^2 + 2*d2(2)*qd2(1)*cos(q2(2))*qd2(2) ...
        + 2*d2(2)*qdd2(1)*sin(q2(2)) + d2(2)*qdd2(2)*sin(q2(2))) ...
        + qdd2(1)*zz21 + qdd2(1)*zz22 + qdd2(2)*zz22 + m22*qdd2(1)*d2(2)^2;
    
    Tt2(2) = (my22*sin(q2(1) + q2(2)) - mx22*cos(q2(1) + q2(2)))*g ...
        + zz22*(qdd2(1) + qdd2(2)) + fv2(2)*qd2(2) + fs2(2)*sign(qd2(2)) ...
        + mx22*(d2(2)*sin(q2(2))*qd2(1)^2 + d2(2)*qdd2(1)*cos(q2(2))) ...
        - my22*(-d2(2)*cos(q2(2))*qd2(1)^2 + d2(2)*qdd2(1)*sin(q2(2)));
   
    % End effector like a punctual mass fixed on B22
    w0p = m4 * [xdd; ydd-g; zeros(4,1)];
    w0r = Psit' * w0p;
    % IDM
    Tta = [Tt1(1); Tt2(1)];
    Ttd = [Tt1(2:3)'; Tt2(2)];

    Jd = (Jtd) \ (Jt * J - Jta);

    T = Tta + J' * w0r + Jd' * Ttd;
