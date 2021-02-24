% TODO: seleccionar modos de ensablaje
function [p] = fkine5(q)
% This function implement the forward kinematic model of a five bar robot
% @parm q: joint positions vector 5x1
% @return p: p end-effector coordinate vector 2X1
    l = 0.205;
    d1 = [-0.125, l, l];
    d2 = [0.125, l, l];
    q1 = q(1:3); % [0.4072, 0, 0];
    q2 = q(4:5); % [1.1393, 0, 0];

    A12 = [d1(2)*cos(q1(1))+d1(1); d1(2)*sin(q1(1))];
    A22 = [d2(2)*cos(q2(1))+d2(1); d2(2)*sin(q2(1))];
    a1 = - 2*d1(1) - 2*d1(2)*cos(q1(1));
    a2 = - 2*d2(1) - 2*d2(2)*cos(q2(1));
    b1 = - 2*sin(q1(1))*d1(2);
    b2 = - 2*sin(q2(1))*d2(2);
    c1 = d1(1)^2 + 2*cos(q1(1))*d1(1)*d1(2) + d1(2)^2 - d1(3)^2;
    c2 = d2(1)^2 + 2*cos(q2(1))*d2(1)*d2(2) + d2(2)^2 - d2(3)^2;

    if b1 ~= b2
        a = a2 - a1;
        b = b2 - b1;
        c = c2 - c1;
        e1 = - a / b;
        e2 = - c / b;
        f1 = 1 + e1^2;
        f2 = 2 * e1 * e2 + a1 + b1 * e1;
        f3 = e2^2 + b1 * e2 + c1;

        x = (-f2 - sqrt(f2^2 - 4 * f1 * f3)) / (2 * f1);
        y = e1 * x + e2;

        A13=[A12(1)-x;A12(2)-y];
        A23=[A22(1)-x;A22(2)-y];
        if  A13(1)*A23(2)-A23(1)*A13(2) <= 0
            x = (-f2 + sqrt(f2^2 - 4 * f1 * f3)) / (2 * f1);
            y = e1 * x + e2;
        end
    else
        x = - (c2 - c1) / (a2 - a1);
        y = (-b1 + (b1^2 - 4 * (x^2 + a1 * x + c1))^(1/2)) / 2;

        A13=[A12(1)-x;A12(2)-y];
        A23=[A22(1)-x;A22(2)-y];
        if  A13(1)*A23(2)-A23(1)*A13(2) <= 0
            y = (-b1 - (b1^2 - 4 * (x^2 + a1 * x + c1))^(1/2)) / 2;
        end
    end

    % for q12 and q22 from a12/a11 and a22/a21
    q1(2) = atan2(y-d1(2)*sin(q1(1)), x-d1(1)-d1(2)*cos(q1(1))) - q1(1);
    q2(2) = atan2(y-d2(2)*sin(q2(1)), x-d2(1)-d2(2)*cos(q2(1))) - q2(1);
    q1(3) = q2(1) + q2(2) - q1(1) - q1(2);
    p = [x;y];

