function [t_sol, u_sol] = ode_mario(metodo, funcion, rango, u0, h)
% Función que resuelve ODEs por el método de runge-kutta
% 
% Parámetros:
%   -método:
%         - 'euler'     Euler method
%         - 'euler+'    Inproved euler method / Heun method
%         - 'b-euler'   Backward euler method
%         - 'midpoint'  Midpoint method 
%         - 'rk3'       Runge-Kutta method 3º order
%         - 'rk4'       Runge-Kutta method 3º order
%         - 'rkf'       Runge-Kutta-Fehlberg method
%   -funcion: Se introduce la EDO de la siguiente forma
%       function dy = funcion(u, t)
%           du = t + 1/u;
%       end
%   -rango: vector que indica el rango de tiempos
%   -x0: valor inicial
%   -h: paso
% 
% Devuelve:
%     - t_sol: vector de tiempos, donde se ha evaluado la función
%     - u_sol: vector de soluciones aproximadas por el método numérico elegido

if nargin < 3
    rango = [0, 4];
    u0 = 0;
    h = 0.1;
elseif nargin < 4
    u0 = 0;
    h = 0.1;
elseif nargin < 5
    h = 0.1;
elseif nargin > 5
    error('Demasiados parámetros')
end

if rango(1) < rango(2)
    t_i = rango(1);
    t_f = rango(2);
elseif rango(1) > rango(2)
    t_i = rango(1);
    t_f = rango(2);
else
    error('Los valores del rango son iguales')
end

t = t_i :h: t_f;

u = zeros(1,length(t));

u(1) = u0;

switch metodo
    case 'euler'
        for contador = 1:length(t)-1        
            u(contador + 1) = u(contador) + h * funcion(u(contador), t(contador));
        end
    case 'euler+'
        for contador = 1:length(t)-1
             k1 = h * funcion(u(contador), t(contador));
             k2 = h * funcion(u(contador) + k1, t(contador + 1));
             u(contador + 1) = u(contador) + 0.5*(k1 + k2);
        end
    case 'b-euler'
%       u_{k+1} == u_k + h * f(u_{k+1}, t_{k+1})
%       0 == u_k + h * f(U, t_{k+1}) - U
        for contador = 1:length(t)-1
            u(contador + 1) = fzero(@(U) u(contador) + funcion(U, t(contador+1)) ...
                * h - U, u(contador));            
        end 
        warning('Poco estable')
    case 'midpoint'
        for contador = 1:length(t)-1
             k1 = h * funcion(u(contador), t(contador));
             k2 = h * funcion(u(contador) + k1 / 2 + k1, t(contador) + h/2);
             u(contador + 1) = u(contador) + k2;
        end         
    case 'rk3'
        for contador = 1:length(t)-1
             k1 = h * funcion(u(contador), t(contador));
             k2 = h * funcion(u(contador) + k1/2, t(contador) + h/2);
             k3 = h * funcion(u(contador) - k1 + 2 * k2, t(contador + 1));
             u(contador + 1) = u(contador) + (1/6) * (k1 + 4*k2 + k3);
        end
    case 'rk4'
        for contador = 1:length(t)-1
             k1 = h * funcion(u(contador), t(contador));
             k2 = h * funcion(u(contador) + k1/2, t(contador) + h/2);
             k3 = h * funcion(u(contador) + k2/2, t(contador) + h/2);
             k4 = h * funcion(u(contador) + k3, t(contador + 1));
             u(contador + 1) = u(contador) + (1/6) * (k1 + 2*k2 +2*k3 + k4);  
        end
    case 'rkf'
        for contador = 1:length(t)-1
             k1 = h * funcion(u(contador), t(contador));
             k2 = h * funcion(u(contador) + k1 / 4, t(contador) + h/4);
             k3 = h * funcion(u(contador) + 3* k1 / 32 + 9* k2/ 32, t(contador) + ...
                 3* h/ 8);
             k4 = h * funcion(u(contador) + 1932* k1/2197 - 7200 * k2 /2197 ...
                 + 7296 * k3 /2197,  t(contador)+ 12*h/13); 
             k5 = h * funcion(u(contador) + 439 * k1 / 216 -8*k2 + 3680 * k3 /513 ...
                 - 845 * k4 / 4104, t(contador) + h);
             k6 = h * funcion(u(contador) - 8 *k1 / 27 + 2*k2 -3544 * k3 /2565 ...
                 + 1859 * k4 / 4104 - 11*k5/40, t(contador) +0.5*h);
             u(contador + 1) = u(contador) + 16/135 * k1 + 0* k2 + 6656/12825 * k3 ...
                 + 28561/56430 *k4 - 9/50 * k5 + 2/55 * k6;
        end        
    otherwise
        warning('Método introducido no válido');
        error('Revisa si está bien escrito');
end
% Parametros que devuelve la función
t_sol = t;
u_sol = u;

% Representación gráfica
figure;
plot(t_sol,u_sol,'Color', [0.4, 0.6, 0.2] ,'LineWidth',1)
title(['Método de ' metodo ' \Deltah= ' num2str(h)])
xlabel('t')
ylabel('u(t)')
grid on

end