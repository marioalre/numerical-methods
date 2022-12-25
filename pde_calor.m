function [vt,vsol] = pde_calor(alpha, xrange, yrange, trange, nx, ny, nt)
%pde_calor Ecuación de calor 2D: ut = alpha(uxx + uyy)
%   Para cambiar las CI y CF, modificar el código
%   [vt, vsol] = pde_calor(1.2,[0,100], [0, 100], [0, 40], 20, 20, 400)
close all
% Condición inicial
A0 = zeros(nx, ny, nt);

val = 0.2;
for ii = 2:nx-1
    m = length(A0(ii:end-ii, 1, 1));
    n = length(A0(1, ii:end-ii, 1));
    A0(ii:end-ii, ii:end-ii, 1) = val*ones(m, n);
    val = val + 0.2;
end

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
for tt = 2:nt
    for ii = 2:nx-1
        for jj = 2:ny-1
            vsol(ii, jj, tt) = vsol(ii,jj,tt-1)*(1-2*rx-2*ry) + ...
                vsol(ii+1, jj, tt-1)*rx + vsol(ii-1, jj, tt-1)*rx + ...
                vsol(ii, jj+1, tt-1)*ry + vsol(ii, jj-1, tt-1)*ry;
        end
    end

    % Para las condiciones iniciales
    % Suponemos condiciones Dirichlet

%     for ii = 1:nx
%         if ii == 1 || ii == nx
%             for jj = 1:ny
%                 if jj == 1 || jj == ny
%                     vsol(ii, jj, tt) = vsol(ii, jj,  tt-1);
%                 end
%             end
%         end
%     end

    % Neumman
%     for ii = 1:nx
%         for jj = 1:ny
%             if ii == 1
%                 vsol(ii, jj, tt) = vsol(ii+1, jj,  tt-1);
%             elseif  ii == nx
%                 vsol(ii, jj, tt) = vsol(ii-1, jj,  tt-1);
%             elseif jj == 1
%                 vsol(ii, jj, tt) = vsol(ii, jj+1,  tt-1);
%             elseif jj == ny
%                 vsol(ii, jj, tt) = vsol(ii, jj-1,  tt-1);
%             end
%         end
%     end

% Mezcla de condiciones
    for ii = 1:nx
        for jj = 1:ny
            if ii == 1
                vsol(ii, jj, tt) = 0;
            elseif  ii == nx
                vsol(ii, jj, tt) = 0;
            elseif jj == 1
                vsol(ii, jj, tt) = vsol(ii, jj+1, tt-1);
            elseif jj == ny
                vsol(ii, jj, tt) = 0;
            end
        end
    end
end
    vt = linspace(trange(1), trange(2), nt);
    
    % Representación gráfica
    figure
    for tt = 1:nt
        surf(vsol(:, : , tt))
        title('Time: '+string(vt(tt))+' s')
        zlim([0, 1.8])
        pause(1e-4)
    end

end