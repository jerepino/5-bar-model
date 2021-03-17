function [t0rd] = fsecondkine5(qdda,qd,q,Ar,J)
% This function implement the second-order forward kinematics model of a five bar robot
% @param qdda: active acceleration joints vector 2X1
% @param qd: velocity joints vector 5x1
% @param q: joint positions vector 5x1
% @param Ar: matrix 2x2 for computation of the dynamic model of the PKM
% @param J: kinematic Jacobian matrix of the PKM
% @return t0rd: derivate twist end-effector vector 2X1
    l = 0.205;
    d1 = [-0.125, l, l];
    d2 = [0.125, l, l];
    q1 = q(1:3);
    q2 = q(4:5);
    qd1 = qd(1:3);
    qd2 = qd(4:5);

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

    at = Ar\b0p;
    t0rd = (J * qdda + at);

