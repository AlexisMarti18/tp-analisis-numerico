% =========================================================================
% ANÁLISIS NUMÉRICO - TRABAJOS DE LABORATORIO
% EJERCICIO 2: Filtrado de una señal muestreada por convolución con un
%              núcleo sinc (filtro pasa bajos ideal truncado)
%
% Uso:  ej2                      -> valores por defecto
%       ej2(dt, Nh)              -> cambia el paso de muestreo y el núcleo
%       ej2(dt, Nh, Fc, Nx)      -> además la frecuencia de corte y Nx
%       [ny, y] = ej2(...)       -> devuelve la respuesta
% =========================================================================
function [ny, y, nx, x, nh, h] = ej2(dt, Nh, Fc, Nx)

  if nargin < 1 || isempty(dt), dt = 0.0001; end  % paso de muestreo [s]
  if nargin < 2 || isempty(Nh), Nh = 101;    end  % largo del núcleo (impar)
  if nargin < 3 || isempty(Fc), Fc = 500;    end  % corte [Hz], 150 < Fc < 1000
  if nargin < 4 || isempty(Nx), Nx = 401;    end  % muestras de x[n]

  if mod(Nh, 2) == 0
    error('Nh debe ser impar para que h[n] quede centrado en n = 0');
  end

  % --- Entrada: muestras de x(t) = 10cos(2pi150t) + 5sen(2pi1000t) -------
  nx = 0:Nx-1;                 % inicio en t = 0
  t  = nx * dt;
  x  = 10*cos(2*pi*150*t) + 5*sin(2*pi*1000*t);

  % --- Respuesta al impulso: h[n] = 2 Fc dt sinc(2 Fc dt n), centrada -----
  M  = (Nh - 1) / 2;
  nh = -M:M;                   % de -(Nh-1)/2 a (Nh-1)/2
  h  = 2*Fc*dt * sinc(2*Fc*dt*nh);

  % --- Respuesta del sistema ----------------------------------------------
  y  = conv(x, h);
  ny = (nx(1) + nh(1)) : (nx(end) + nh(end));

  % --- Gráficos: misma escala, límites = elementos no nulos de y[n] -------
  lim = [ny(1), ny(end)];

  figure('Name', sprintf('Ejercicio 2 - dt = %g s, Nh = %d, Fc = %g Hz', dt, Nh, Fc));

  subplot(3, 1, 1);
  stem(nx, x, 'filled', 'markersize', 2);
  xlim(lim); grid on;
  title(sprintf('Entrada x[n]  (N_x = %d, \\Deltat = %g s, F_m = %g Hz)', Nx, dt, 1/dt));
  ylabel('x[n]');

  subplot(3, 1, 2);
  stem(nh, h, 'filled', 'markersize', 2);
  xlim(lim); grid on;
  title(sprintf('Respuesta al impulso h[n] = 2F_c\\Deltat sinc(2F_c\\Deltat n)  (N_h = %d, F_c = %g Hz)', Nh, Fc));
  ylabel('h[n]');

  subplot(3, 1, 3);
  stem(ny, y, 'filled', 'markersize', 2);
  xlim(lim); grid on;
  title('Respuesta y[n] = x[n] * h[n]');
  xlabel('n'); ylabel('y[n]');

end
