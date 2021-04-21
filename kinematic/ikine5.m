% TODO: Selecionar working mode
function [q,qa,qp] = ikine5(p, q_prev)
% This function implement the inverse kinematic model of a five bar robot
% @parm p: p end-effector coordinate vector 2X1
% @parm q_prev: previos joint position vector 5X1
% @return q: joint positions vector 5x1
% @return qa: active joints positions vector 2X1
% @return qp: underactuated joints positions vector 3X1

    l = 0.205;
    d11 = -0.125;
    d12 = l;
    d13 = l;
    d21 = 0.125;
    d22 = l;
    d23 = l;
    x = p(1);
    y = p(2);
    q_aux = zeros(5,4);

    dist013 = sqrt((p(1)-d11)^2+p(2)^2);
    dist023 = sqrt((p(1)-d21)^2+p(2)^2);
    if dist013 <= 2*l && dist023 <= 2*l
        % q11 = -2*atan((((d11^2 - 2*d11*x - d12^2 + 2*d12*d13 - d13^2 + x^2 + y^2)*(- d11^2 + 2*d11*x + d12^2 + 2*d12*d13 + d13^2 - x^2 - y^2))^(1/2) - 2*d12*y)/(d11^2 - 2*d11*d12 - 2*d11*x + d12^2 + 2*d12*x - d13^2 + x^2 + y^2));
        a1 = -2*d12*(x-d11);
        b1 = -2*d12*y;
        c1 = (x-d11)^2+y^2+d12^2-d13^2;
        q11 = 2*atan((-b1-sqrt(b1^2-c1^2+a1^2)) / (c1-a1));
        if ~isreal(q11)
            q_aux(:,1) = [NaN, NaN, NaN, NaN, NaN]';
            q_aux(:,2) = [NaN, NaN, NaN, NaN, NaN]';
        else
            % q21 = -2*atan((((d21^2 - 2*d21*x - d22^2 + 2*d22*d23 - d23^2 + x^2 + y^2)*(- d21^2 + 2*d21*x + d22^2 + 2*d22*d23 + d23^2 - x^2 - y^2))^(1/2) - 2*d22*y)/(d21^2 - 2*d21*d22 - 2*d21*x + d22^2 + 2*d22*x - d23^2 + x^2 + y^2));
            a2 = -2*d22*(x-d21);
            b2 = -2*d22*y;
            c2 = (x-d21)^2+y^2+d22^2-d23^2;
            q21 = 2*atan((-b2+sqrt(b2^2-c2^2+a2^2)) /(c2-a2));
            if ~isreal(q21)
                q_aux(:,1) = [NaN, NaN, NaN, NaN, NaN]';
            else
                q12 = - q11 + atan2(y - d12*sin(q11), x - d11 - d12*cos(q11));
                q22 = - q21 + atan2(y - d22*sin(q21), x - d21 - d22*cos(q21));
                q13 = q21 - q12 - q11 + q22;
                q_aux(:,1) = [q11, q12, q13, q21, q22]';
            end
            % q21 = 2*atan((((d21^2 - 2*d21*x - d22^2 + 2*d22*d23 - d23^2 + x^2 + y^2)*(- d21^2 + 2*d21*x + d22^2 + 2*d22*d23 + d23^2 - x^2 - y^2))^(1/2) + 2*d22*y)/(d21^2 - 2*d21*d22 - 2*d21*x + d22^2 + 2*d22*x - d23^2 + x^2 + y^2));
            q21 = 2*atan((-b2-sqrt(b2^2-c2^2+a2^2)) /(c2-a2));
            if ~isreal(q21)
                q_aux(:,2) = [NaN, NaN, NaN, NaN, NaN]';
            else
                q12 = - q11 + atan2(y - d12*sin(q11), x - d11 - d12*cos(q11));
                q22 = - q21 + atan2(y - d22*sin(q21), x - d21 - d22*cos(q21));
                q13 = q21 - q12 - q11 + q22;
                q_aux(:,2) = [q11, q12, q13, q21, q22]';
            end
        end

        % q11 = 2*atan((((d11^2 - 2*d11*x - d12^2 + 2*d12*d13 - d13^2 + x^2 + y^2)*(- d11^2 + 2*d11*x + d12^2 + 2*d12*d13 + d13^2 - x^2 - y^2))^(1/2) + 2*d12*y)/(d11^2 - 2*d11*d12 - 2*d11*x + d12^2 + 2*d12*x - d13^2 + x^2 + y^2));
        q11 = 2*atan((-b1+sqrt(b1^2-c1^2+a1^2)) / (c1-a1));
        if ~isreal(q11)
            q_aux(:,3) = [NaN, NaN, NaN, NaN, NaN]';
            q_aux(:,4) = [NaN, NaN, NaN, NaN, NaN]';
        else
            % q21 = -2*atan((((d21^2 - 2*d21*x - d22^2 + 2*d22*d23 - d23^2 + x^2 + y^2)*(- d21^2 + 2*d21*x + d22^2 + 2*d22*d23 + d23^2 - x^2 - y^2))^(1/2) - 2*d22*y)/(d21^2 - 2*d21*d22 - 2*d21*x + d22^2 + 2*d22*x - d23^2 + x^2 + y^2));
            q21 = 2*atan((-b2+sqrt(b2^2-c2^2+a2^2)) /(c2-a2));
            if ~isreal(q21)
                q_aux(:,3) = [NaN, NaN, NaN, NaN, NaN]';
            else
                q12 = - q11 + atan2(y - d12*sin(q11), x - d11 - d12*cos(q11));
                q22 = - q21 + atan2(y - d22*sin(q21), x - d21 - d22*cos(q21));
                q13 = q21 - q12 - q11 + q22;
                q_aux(:,3) = [q11, q12, q13, q21, q22]';
            end
            % q21 = 2*atan((((d21^2 - 2*d21*x - d22^2 + 2*d22*d23 - d23^2 + x^2 + y^2)*(- d21^2 + 2*d21*x + d22^2 + 2*d22*d23 + d23^2 - x^2 - y^2))^(1/2) + 2*d22*y)/(d21^2 - 2*d21*d22 - 2*d21*x + d22^2 + 2*d22*x - d23^2 + x^2 + y^2));
            q21 = 2*atan((-b2-sqrt(b2^2-c2^2+a2^2)) /(c2-a2));
            if ~isreal(q21)
                q_aux(:,4) = [NaN, NaN, NaN, NaN, NaN]';
            else
                q12 = - q11 + atan2(y - d12*sin(q11), x - d11 - d12*cos(q11));
                q22 = - q21 + atan2(y - d22*sin(q21), x - d21 - d22*cos(q21));
                q13 = q21 - q12 - q11 + q22;
                q_aux(:,4) = [q11, q12, q13, q21, q22]';
            end
        end
        % find the correct working mode base in the previous mode
        % qa(:,i) = [q(1,3);q(4,3)]; 
        if nargin < 2
            q = q_aux; %(:,1);
            qa = [q(1,:); q(4,:)];
            qp = [q(2:3,:); q(5,:)];
        else
            distance_between_q = sum(angdiff(q_aux, q_prev * ones(1,4)).^2);
            index = find(distance_between_q == min(distance_between_q));
            q = q_aux(:,index(1));
            qa = [q(1); q(4)];
            qp = [q(2:3); q(5)];
        end
    else
        disp('ERROR: point out of work space');
            q = [NaN; NaN; NaN; NaN; NaN];
            qa = [q(1); q(4)];
            qp = [q(2:3); q(5)];
        % point out of work space
    end

