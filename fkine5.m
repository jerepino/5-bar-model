function [p,q] = fkine5(qa, assembly)
% This function implement the forward kinematic model of a five bar robot
% @parm qa: active joint positions vector 2x1
% @parm assembly: assembly mode selecction 1 for + -1 for -
% @return p: p end-effector coordinate vector 2X1
% @return q: joint positon vector 5X1
    l = 0.205;
    d1 = [-0.125, l, l];
    d2 = [0.125, l, l];
    q11 = qa(1);
    q21 = qa(2);

    A12 = [d1(2)*cos(q11)+d1(1); d1(2)*sin(q11)];
    A22 = [d2(2)*cos(q21)+d2(1); d2(2)*sin(q21)];

    if sqrt((A12-A22)' * (A12-A22)) <= 2*l
        a1 = - 2*d1(1) - 2*d1(2)*cos(q11);
        a2 = - 2*d2(1) - 2*d2(2)*cos(q21);
        b1 = - 2*sin(q11)*d1(2);
        b2 = - 2*sin(q21)*d2(2);
        c1 = d1(1)^2 + 2*cos(q11)*d1(1)*d1(2) + d1(2)^2 - d1(3)^2;
        c2 = d2(1)^2 + 2*cos(q21)*d2(1)*d2(2) + d2(2)^2 - d2(3)^2;

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

            A13=[x-A12(1);y-A12(2)];
            A23=[x-A22(1);y-A22(2)];
            if  A13(1)*A23(2)-A23(1)*A13(2) <= 0 && assembly == 1
                x = (-f2 + sqrt(f2^2 - 4 * f1 * f3)) / (2 * f1);
                y = e1 * x + e2;
            elseif  A13(1)*A23(2)-A23(1)*A13(2) >= 0 && assembly == -1
                x = (-f2 + sqrt(f2^2 - 4 * f1 * f3)) / (2 * f1);
                y = e1 * x + e2;
            end
        else
            x = - (c2 - c1) / (a2 - a1);
            y = (-b1 + (b1^2 - 4 * (x^2 + a1 * x + c1))^(1/2)) / 2;

            A13=[x-A12(1);y-A12(2)];
            A23=[x-A22(1);y-A22(2)];
            if  A13(1)*A23(2)-A23(1)*A13(2) <= 0 && assembly == 1
                y = (-b1 - (b1^2 - 4 * (x^2 + a1 * x + c1))^(1/2)) / 2;
            elseif  A13(1)*A23(2)-A23(1)*A13(2) >= 0 && assembly == -1
                y = (-b1 - (b1^2 - 4 * (x^2 + a1 * x + c1))^(1/2)) / 2;
            end
        end

        % for q12 and q22 from a12/a11 and a22/a21
        q12 = atan2(y-d1(2)*sin(q11), x-d1(1)-d1(2)*cos(q11)) - q11;
        q22 = atan2(y-d2(2)*sin(q21), x-d2(1)-d2(2)*cos(q21)) - q21;
        q13 = q21 + q22 - q11 - q12;
        p = [x;y];
        q = [q11; q12; q13; q21; q22];

    else
        disp('ERROR: Euclidian distance > 2*l');
        p = [NaN; NaN];
        q = [NaN; NaN; NaN; NaN; NaN];
    end
