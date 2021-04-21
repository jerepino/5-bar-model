function [qd,qda,qdu,Ar,B,J,Jinv,Jt,Jta,Jtd,Psit] = ifirstkine5(t0r,p,q)
% This function implement the firs-order inverse kinematics model of a five bar robot
% @param t0r: twist end-effector vector 2X1
% @param p: end-effector coordinate vector 2X1
% @param q: joint positions vector 5x1
% @return qd: velocity joints vector 5x1
% @return qda: active velocity joints vector 2X1
% @return qdu: underactuated velocity joints vector 3X1
% @return Ar: matrix 2x2 relates the n dof independent coordinates of the platform twist 0 t r to the actuated joints q̇ a
% @return B: matrix 2x2 for computation of the dynamic model of the PKM
% @return J: kinematic Jacobian matrix of the PKM
% @return Jinv: inverse kinematic Jacobian matrix of the PKM
% @return Jt: crucial for computation of the dynamic model of the PKM.
% @return Jta: crucial for computation of the dynamic model of the PKM.
% @return Jtd: crucial for computation of the dynamic model of the PKM.
% @return Psit: matrix 2X6

    l = 0.205;
    d1 = [-0.125, l, l];
    d2 = [0.125, l, l];
    q1 = q(1:3);
    q2 = q(4:5);
    Phi = sum(q1);
    x = p(1);
    y = p(2);

    B = -[ d1(2)*sin(q1(2)), 0;
          0, d2(2)*sin(q2(2))];


    Ar = [cos(q1(1) + q1(2)), sin(q1(1) + q1(2));
          cos(q2(1) + q2(2)), sin(q2(1) + q2(2))];
    % forward model
    J = -Ar \ B ;
    % t0r = J * qda;
    % t0r =[ -(d1(2)*qd11*sin(q2(1) + q2(2))*sin(q1(2)) - ...
    %                 d2(2)*qd21*sin(q1(1) + q1(2))*sin(q2(2)))/sin(q1(1) + q1(2) - q2(1) - q2(2));
    %         (d1(2)*qd11*cos(q2(1) + q2(2))*sin(q1(2)) - ...
    %                 d2(2)*qd21*cos(q1(1) + q1(2))*sin(q2(2)))/sin(q1(1) + q1(2) - q2(1) - q2(2))];

    Jinv = -B \ Ar;
    qda = Jinv * t0r;
    % qda = [(xd*cos(q1(1) + q1(2)) + yd*sin(q1(1) + q1(2)))/(d1(2)*sin(q1(2)));
    %        (xd*cos(q2(1) + q2(2)) + yd*sin(q2(1) + q2(2)))/(d2(2)*sin(q2(2)))];

    % passive velocities calculation

    ay = cos(Phi)^2 / (x - d2(2)*cos(q2(1)) - d2(1));
    ax = - ((y-d2(2)*sin(q2(1)))*cos(Phi)^2) ...
           / (x - d2(2)*cos(q2(1)) - d2(1))^2;
    aq = -(d2(2) * ((x-d2(1))*cos(q2(1)) ...
          +y*sin(q2(1))-d2(2))*cos(Phi)^2) ...
          / (x - d2(2)*cos(q2(1)) - d2(1))^2;

    Psit = [1                     0
            0                     1
            0                     0
            0                     0
            0                     0
            (ax+aq*Jinv(2,1)) (ay+aq*Jinv(2,2))];

    Jt = [-sin(q1(1) + q1(2)),cos(q1(1) + q1(2));
            (ax + aq*Jinv(2,1)), (ay + aq*Jinv(2,2));
            -sin(q2(1) + q2(2)),cos(q2(1) + q2(2))];


    Jta = [ d1(3) + d1(2)*cos(q1(2)), 0;
            1,                0;
            0,      d2(3) + d2(2)*cos(q2(2))];

    Jtd = [ d1(3), 0,   0;
            1,   1,   0;
            0,   0, d2(3)];

    qdu = Jtd\(Jt * t0r - Jta * qda);
    qd = [qda(1);qdu(1);qdu(2);qda(2);qdu(3)];
    % qdu = [ - (xd*cos(q1(1) + q1(2)) + yd*sin(q1(1) + q1(2)))/(d1(2)*sin(q1(2))) - ...
    %                 (d1(2)*xd*cos(q1(1)) + d1(2)*yd*sin(q1(1)))/(d1(2)*d1(3)*sin(q1(2)));
    %          (d2(2)*xd*cos(q1(1)) + d2(2)*yd*sin(q1(1)) + ax*d1(3)*d2(2)*xd*sin(q1(2)) + ...
    %                 ay*d1(3)*d2(2)*yd*sin(q1(2)))/(d1(3)*d2(2)*sin(q1(2))) + ...
    %                 (aq*d1(3)*xd*cos(q2(1) + q2(2))*sin(q1(2)) + ...
    %                 aq*d1(3)*yd*sin(q2(1) + q2(2))*sin(q1(2)))/(d1(3)*d2(2)*sin(q1(2))*sin(q2(2)));
    %         - (xd*cos(q2(1) + q2(2)) + yd*sin(q2(1) + q2(2)))/(d2(2)*sin(q2(2))) - ...
    %                 (d2(2)*xd*cos(q2(1)) + d2(2)*yd*sin(q2(1)))/(d2(2)*d2(3)*sin(q2(2)))];
