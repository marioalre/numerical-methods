function [vt,vsol] = pde_ondas(alpha, xrange, yrange, trange, nx, ny, nt)
%pde_ondas Ecuación de ondas 2D: utt = alpha(uxx + uyy)
%   Para cambiar las CI y CF, modificar el código
%   [vt, vsol] = pde_ondas(0.8,[0,200], [0, 200], [0, 10], 100, 100, 1000);
close all

% Condición inicial
A0 = zeros(nx, ny, nt+1);

val = 1;
for ii = 1:nx-1
    if ii > 20
        m = length(A0(ii:end-ii, 1, 2)); % Dejamos el primer elemento para representar las velocidades
        n = length(A0(1, ii:end-ii, 2));
        A0(ii:end-ii, ii:end-ii, 2) = val*ones(m, n);
        val = val + 0.6;
    end
end

% figure
% surf(A0)
% title('Condición inicial')

% CF Neuman
cf_left = 0;
cf_right = 0;
cf_top = 0;
cf_bottom = 0;

% CF Dirichlet, sumar a todos lo elementos de la línea 8 el valor deseado

% Inicializando programa
vsol = A0;
xvector = linspace(xrange(1), xrange(2), nx);
dx = xvector(2) - xvector(1);

yvector = linspace(yrange(1), yrange(2), nx);
dy = yvector(2) - yvector(1);

tvector = linspace(trange(1), trange(2), nx);
dt = tvector(2) - tvector(1);

rx = alpha*dt/dx^2;
ry = alpha*dt/dy^2;

if  rx>0.5 || ry>0.5
    warning('Método inestable')
    disp(rx)
    disp(ry)
end

% Stencil
for tt = 3:nt
    for ii = 2:nx-1
        for jj = 2:ny-1
            if tt == 3
            vsol(ii, jj, tt) = vsol(ii,jj,tt-1)*(1-rx-ry) + ...
                0.5*vsol(ii+1, jj, tt-1)*rx + 0.5*vsol(ii-1, jj, tt-1)*rx + ...
                0.5*vsol(ii, jj+1, tt-1)*ry + 0.5*vsol(ii, jj-1, tt-1)*ry + ...
                dt*vsol(ii, jj, tt-2);                
            else
            vsol(ii, jj, tt) = vsol(ii,jj,tt-1)*(2-2*rx-2*ry) + ...
                vsol(ii+1, jj, tt-1)*rx + vsol(ii-1, jj, tt-1)*rx + ...
                vsol(ii, jj+1, tt-1)*ry + vsol(ii, jj-1, tt-1)*ry - ...
                vsol(ii, jj, tt-2);
            end
        end
    end

    % Para las condiciones iniciales
    % Suponemos condiciones Dirichlet

%     for ii = 2:nx
%         if ii == 1 || ii == nx
%             for jj = 2:ny
%                 if jj == 1 || jj == ny
%                     vsol(ii, jj, tt) = vsol(ii, jj,  tt-1);
%                 end
%             end
%         end
%     end

    % Neumman
    for ii = 1:nx
        for jj = 1:ny
            if ii == 1
                vsol(ii, jj, tt) = vsol(ii+1, jj,  tt-1) + dt*cf_left;
            elseif  ii == nx
                vsol(ii, jj, tt) = vsol(ii-1, jj,  tt-1) + dt*cf_right;
            elseif jj == 1
                vsol(ii, jj, tt) = vsol(ii, jj+1,  tt-1) + dt*cf_bottom;
            elseif jj == ny
                vsol(ii, jj, tt) = vsol(ii, jj-1,  tt-1) + dt*cf_top;
            end
        end
    end


end
    vt = linspace(trange(1), trange(2), nt);
    
    % Representación gráfica animada
    figure
    for tt = 1:nt
        surf(vsol(:, : , tt))
        title('Time: '+string(vt(tt))+' s')
        zlim([-8, 12])
%         colorbar()
        pause(1e-2)
    end

end