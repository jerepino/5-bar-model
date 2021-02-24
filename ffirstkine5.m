function [t0r] = ffirstkine5(q,qda)
% This function implement the firs-order forward kinematics model of a five bar robot
% @parm q: joint positions vector 5x1
% @parm qda: active velocity joints vector 2X1
% @return t0r: twist end-effector vector 2X1
    l = 0.205;
    d1 = [-0.125, l, l];
    d2 = [0.125, l, l];
    q1 = q(1:3);
    q2 = q(4:5);

    B = [ -d1(2)*sin(q1(2)), 0;
          0, -d2(2)*sin(q2(2))];

    Ar = [cos(q1(1) + q1(2)), sin(q1(1) + q1(2));
          cos(q2(1) + q2(2)), sin(q2(1) + q2(2))];
    % forward model
    J = -B / Ar;
    t0r = J * qda;
    % qd11 = qda(1);
    % qd21 = qda(2);
    % t0r =[ -(d1(2)*qd11*sin(q2(1) + q2(2))*sin(q1(2)) - ...
    %                 d2(2)*qd21*sin(q1(1) + q1(2))*sin(q2(2)))/sin(q1(1) + q1(2) - q2(1) - q2(2));
    %         (d1(2)*qd11*cos(q2(1) + q2(2))*sin(q1(2)) - ...
    %                 d2(2)*qd21*cos(q1(1) + q1(2))*sin(q2(2)))/sin(q1(1) + q1(2) - q2(1) - q2(2))];


