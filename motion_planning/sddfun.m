function [sdd] = sddfun(t,T,tau)
%%
%   Trapezoidal timing law
%   Funcion es la derivada de ds  y son 2 escalones
%   Esta funcion se puede suavizar para tener jerk constante!!
%%
    a = 1 / (T * tau);
    if t < 0
        sdd = 0;
    elseif t <= tau
        sdd = a;
    elseif t > tau && t <= T
        sdd = 0;
    elseif t > T && t <= (T+tau)
        sdd = - a;
    elseif t > (T+tau)
        sdd = 0;
    end
end