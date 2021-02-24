%% NOTA: Evaluar pasar más parametros como la matris Psit
function [qdd,qdda,qddu,Psit1,Psit2,atp,ad] = isecondkine5(t0rd,t0r,p,q,qd,Ar,B,Jinv,Jt,Jta,Jtd,Psit)
% This function implement the second-order inverse kinematics model of a five bar robot
% @param t0rd: derivate twist end-effector vector 2X1
% @param t0r: twist end-effector vector 2X1
% @param p: end-effector coordinate vector 2X1
% @param q: joint positions vector 5x1
% @param qd: velocity joints vector 5x1
% @param Ar: matrix 2x2 for computation of the dynamic model of the PKM
% @param B: matrix 2x2 for computation of the dynamic model of the PKM
% @param Jinv: inverse kinematic Jacobian matrix of the PKM
% @param Jt: crucial for computation of the dynamic model of the PKM.
% @param Jta: crucial for computation of the dynamic model of the PKM.
% @param Jtd: crucial for computation of the dynamic model of the PKM.
% @return qdd: acceleration joint vector 5x1
% @return qdda: active acceleration joints vector 2X1
% @return qddu: underactuated acceleration joints vector 3X1
% @return Psit1: matrix 2x6
% @return Psit2: matrix 1x6
% @return atp: acceleration using for dynamics
% @return ad: acceleration using for dynamics

    l = 0.205;
    d1 = [-0.125, l, l];
    d2 = [0.125, l, l];
    q1 = q(1:3);
    q2 = q(4:5);
    qd1 = qd(1:3);
    qd2 = qd(4:5);
    Phi = sum(q1);
    Phid = sum(qd1);
    x = p(1);
    y = p(2);
    xd = t0r(1);
    yd = t0r(2);

    wr1 = [cos(q1(1)+q1(2)); sin(q1(1)+q1(2)); 0; 0; 0; 0];
    wr2 = [cos(q2(1)+q2(2)); sin(q2(1)+q2(2)); 0; 0; 0; 0];

    b01a_ = [-d1(2)*qd1(1)*cos(q1(1)) - d1(3)*cos(q1(1)+q1(2))*(qd1(1)+qd1(2));
             -d1(2)*qd1(1)*sin(q1(1)) - d1(3)*sin(q1(1)+q1(2))*(qd1(1)+qd1(2));
             zeros(4,1)];
    b01b_ = [-d1(3)*cos(q1(1)+q1(2));
             -d1(3)*sin(q1(1)+q1(2));
             zeros(4,1)];
    b0p1 = b01a_ * qd1(1) + b01b_ * qd1(2)*(qd1(1)+qd1(2));

    b02a_ = [-d2(2)*qd2(1)*cos(q2(1)) - d2(3)*cos(q2(1)+q2(2))*(qd2(1)+qd2(2));
             -d2(2)*qd2(1)*sin(q2(1)) - d2(3)*sin(q2(1)+q2(2))*(qd2(1)+qd2(2));
             zeros(4,1)];
    b02b_ = [-d2(3)*cos(q2(1)+q2(2));
             -d2(3)*sin(q2(1)+q2(2));
             zeros(4,1)];
    b0p2 = b02a_ * qd2(1) + b02b_ * qd2(2)*(qd2(1)+qd2(2));

    b0p = [(wr1'*b0p1); (wr2'*b0p2)];

    aqq = inv(B) * b0p;
    qdda = (Jinv * t0rd + aqq);

    % Aceleracion pasiva

    ayd = (Phid*sin(2*Phi))/(d2(1) - x + d2(2)*cos(q2(1))) - (cos(Phi)^2*(xd + d2(2)*qd2(1)*sin(q2(1))))/(d2(1) - x + d2(2)*cos(q2(1)))^2;

    axd = (Phid*sin(2*Phi)*(y - d2(2)*sin(q2(1))))/(d2(1) - x + d2(2)*cos(q2(1)))^2 ...
        - (cos(Phi)^2*(yd - d2(2)*qd2(1)*cos(q2(1))))/(d2(1) - x + d2(2)*cos(q2(1)))^2 ...
        - (2*cos(Phi)^2*(xd + d2(2)*qd2(1)*sin(q2(1)))*(y - d2(2)*sin(q2(1))))/(d2(1) - x + d2(2)*cos(q2(1)))^3;

    aqd = (2*d2(2)*cos(Phi)^2*(xd + d2(2)*qd2(1)*sin(q2(1)))*(d2(2) - sin(q2(1))*y + cos(q2(1))*(d2(1) - x)))...
        /(d2(1) - x + d2(2)*cos(q2(1)))^3 ...
        - (d2(2)*cos(Phi)^2*(xd*cos(q2(1)) + yd*sin(q2(1)) + qd2(1)*cos(q2(1))*y + qd2(1)*sin(q2(1))*(d2(1) - x)))...
        /(d2(1) - x + d2(2)*cos(q2(1)))^2 ...
        - (2*Phid*d2(2)*cos(Phi)*sin(Phi)*(d2(2) - sin(q2(1))*y + cos(q2(1))*(d2(1) - x)))...
        /(d2(1) - x + d2(2)*cos(q2(1)))^2;

    j21_invd = -(sin(q2(1) + q2(2))*(qd2(1) + qd2(2)))/(d2(2)*sin(q2(2))) ...
               -(qd2(2)*cos(q2(1) + q2(2))*cos(q2(2)))/(d2(2)*sin(q2(2))^2);

    j22_invd = -(qd2(2)*sin(q2(1) + q2(2))*cos(q2(2)))/(d2(2)*sin(q2(2))^2) ...
               +(cos(q2(1) + q2(2))*(qd2(1) + qd2(2)))/(d2(2)*sin(q2(2)));

    aq = -(d2(2) * ((x-d2(1))*cos(q2(1))+y*sin(q2(1))-d2(2))*cos(Phi)^2) / (x - d2(2)*cos(q2(1)) - d2(1))^2;

    psitd1 = axd + aqd*Jinv(2,1)+ aq*j21_invd;
    psitd2 = ayd + aqd*Jinv(2,2)+ aq*j22_invd;

    Psidt = [0 0; 0 0; 0 0; 0 0; 0 0; psitd1 psitd2];

    Psit1 = [-sin(q1(1)+q1(2)) cos(q1(1)+q1(2)) 0 0 0 0
             0                 0                0 0 0 1];
    Psit2 = [-sin(q2(1)+q2(2)) cos(q2(1)+q2(2)) 0 0 0 0];

    dc1 = Psit1 * (eye(6)* Psidt * t0r - b0p1);
    dc2 = Psit2 * (eye(6)* Psidt * t0r - b0p2);
    dc = [dc1;
          dc2];

    qddu =  inv(Jtd)*(Jt*t0rd - Jta*qdda + dc);

    qdd = [qdda(1);qddu(1);qddu(2);qdda(2);qddu(3)];

    at = inv(Ar) * b0p;
    atp = Psit * at + Psidt * t0r;% using for dynamics
    ad = inv(Jtd) * (dc + Jt*at); % For dynamics
