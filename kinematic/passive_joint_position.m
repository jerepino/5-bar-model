 function [q,qp] = passive_joint_position(qa, p)
 % Function for passive position calculation
    l = 0.205;
    d1 = [-0.125, l, l];
    d2 = [0.125, l, l];
    q11 = qa(1);
    q22 = qa(2);
    x = p(1);
    y = p(2);
    q12 = atan2(y-d1(2)*sin(q11), x-d1(1)-d1(2)*cos(q11)) - q11;
    q22 = atan2(y-d2(2)*sin(q21), x-d2(1)-d2(2)*cos(q21)) - q21;
    q13 = q21 + q22 - q11 - q12;
    q = [q11; q12; q13; q21; q22];
    qp = [q12; q13; q22];

