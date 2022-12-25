% Adams-Moulton Methods
function [vt,vsol]=ODE_AdamsMoulton(t0,h,tfinal,f,u0, tol)

% RK4 formula. u(k+1)=u(k)+h(s1+ 2*s2+ 2*s3 + s4)/6
% u'=f(t,u), t0<t<tfinal, u0=u(t0) Initial Condition
    
    vsol=u0; %  Valor inicial del vector vsol
    u=u0; %     Valor inicial de u
    t = t0; %   Valor inicial del tiempo

    for i=1:3
        s1 = f(t, u);% pendiente para el tiempo t
        s2 = f(t + h/2, u + h*s1/2); 
        s3 = f(t + h/2, u + h*s2 /2);
        s4 = f(t + h, u + h*s3);
        u = u + h/6 *(s1 + 2* s2 + 2* s3 + s4);
        vsol=[vsol; u]; % agregamos al vector solucion
        t = t + h;
    end
    vt=(t0:h:t0+3*h)';

    while vt(end) < tfinal
        y_predict = vsol(end) + h/24*(55*f(vt(end), vsol(end)) - 59*f(vt(end-1), vsol(end-1)) ...
            + 37*f(vt(end-2), vsol(end-2)) - 9*f(vt(end-3), vsol(end-3)));
    
        y_correct = vsol(end) + h/24*(9*f(vt(end) + h, y_predict) + 19*f(vt(end), vsol(end)) ...
            - 5*f(vt(end-1), vsol(end-1)) + f(vt(end-2), vsol(end-2)));
    
        err = 1/15 * (y_correct - y_predict);
        
        % Control del paso
        if abs(err) > tol
            h = h/2;
            disp('Disminuimos el paso')
            t_int = vt(end-2:end);
            v_int = vsol(end-2:end);
            tq_int = vt(end-2):h:vt(end);
            tq_int = tq_int(end-3:end);
            vq_int = interp1(t_int, v_int, tq_int, "spline");

            vsol(end-1:end) = vq_int(end-3:end-2);
            vsol = [vsol; vq_int(end-1); vq_int(end)];

            vt(end-1:end) = tq_int(end-3:end-2);
            vt = [vt; tq_int(end-1); tq_int(end)];
            
            if h > 1e-4
            % Para que rehaga los calculos hasta alcanzar la tolerancia 
                continue
            else
                h = 1e-4; % Para controlar que el paso no vaya a 0
            end        
        end
        fprintf('El valor del paso es de: %f \n', h)

        vsol = [vsol; y_correct];
        vt = [vt; vt(end) + h];
    end
    
    % Representacion grafica
    t = linspace(t0, tfinal, 400);
    u = interp1(vt, vsol, t, 'spline');
    plot(t,u, 'm-', 'LineWidth',1.0);
    hold on
    plot(vt,vsol, 'Color' ,[71, 56, 179]./255, 'Marker','o', 'MarkerFaceColor','auto');
    grid minor
    title('Solución por el método Adams-Moulton')
    xlabel('t')
    ylabel('u(t)')
    xlim([t0, tfinal])
    hold off
end