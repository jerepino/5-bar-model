% TODO: Selecionar working mode
% TODO: verificar si es real, si no devolver algo y avisar
function [q,qa,qu] = ikine5(p)
% This function implement the inverse kinematic model of a five bar robot
% @parm p: p end-effector coordinate vector 2X1
% @return q: joint positions vector 5x1
% @return qa: active joints positions vector 2X1
% @return qu: underactuated joints positions vector 3X1

    l = 0.205;
    d11 = -0.125;
    d12 = l;
    d13 = l;
    d21 = 0.125;
    d22 = l;
    d23 = l;
    x = p(1);
    y = p(2);

    % q11 = -2*atan((((d11^2 - 2*d11*x - d12^2 + 2*d12*d13 - d13^2 + x^2 + y^2)*(- d11^2 + 2*d11*x + d12^2 + 2*d12*d13 + d13^2 - x^2 - y^2))^(1/2) - 2*d12*y)/(d11^2 - 2*d11*d12 - 2*d11*x + d12^2 + 2*d12*x - d13^2 + x^2 + y^2));
    q11 = 2*atan((((d11^2 - 2*d11*x - d12^2 + 2*d12*d13 - d13^2 + x^2 + y^2)*(- d11^2 + 2*d11*x + d12^2 + 2*d12*d13 + d13^2 - x^2 - y^2))^(1/2) + 2*d12*y)/(d11^2 - 2*d11*d12 - 2*d11*x + d12^2 + 2*d12*x - d13^2 + x^2 + y^2));
    q21 = -2*atan((((d21^2 - 2*d21*x - d22^2 + 2*d22*d23 - d23^2 + x^2 + y^2)*(- d21^2 + 2*d21*x + d22^2 + 2*d22*d23 + d23^2 - x^2 - y^2))^(1/2) - 2*d22*y)/(d21^2 - 2*d21*d22 - 2*d21*x + d22^2 + 2*d22*x - d23^2 + x^2 + y^2));
    % q21 = 2*atan((((d21^2 - 2*d21*x - d22^2 + 2*d22*d23 - d23^2 + x^2 + y^2)*(- d21^2 + 2*d21*x + d22^2 + 2*d22*d23 + d23^2 - x^2 - y^2))^(1/2) + 2*d22*y)/(d21^2 - 2*d21*d22 - 2*d21*x + d22^2 + 2*d22*x - d23^2 + x^2 + y^2));
    q12 = - q11 + atan2(y - d12*sin(q11), x - d11 - d12*cos(q11));
    q22 = - q21 + atan2(y - d22*sin(q21), x - d21 - d22*cos(q21));
    q13 = q21 - q12 - q11 + q22;
    q = [q11;q12;q13;q21;q22];
    qa = [q11;q21];
    qu = [q12;q13;q22];
