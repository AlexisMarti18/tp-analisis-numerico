% =========================================================================
% ANÁLISIS NUMÉRICO - TRABAJOS DE LABORATORIO
% EJERCICIO 3: Ejercicio 2 resuelto con la propiedad de convolución de la
%              transformada de Fourier (algoritmo FFT)
%
%   y[n] = x[n] * h[n]   <=>   Y[k] = X[k] . H[k]
%
% Uso:  ej3                     -> valores por defecto
%       ej3(dt, Nh, Fc, Nx, N)
% =========================================================================
function [ny, y_fft, Y] = ej3(dt, Nh, Fc, Nx, N)

  if nargin < 1 || isempty(dt), dt = 0.0001; end
  if nargin < 2 || isempty(Nh), Nh = 101;    end
  if nargin < 3 || isempty(Fc), Fc = 500;    end
  if nargin < 4 || isempty(Nx), Nx = 401;    end
  if nargin < 5 || isempty(N),  N  = 1024;   end   % N = 2^r

  % --- Señales del Ejercicio 2 --------------------------------------------
  nx = 0:Nx-1;
  x  = 10*cos(2*pi*150*nx*dt) + 5*sin(2*pi*1000*nx*dt);

  M  = (Nh - 1) / 2;
  nh = -M:M;
  h  = 2*Fc*dt * sinc(2*Fc*dt*nh);

  % --- Largo necesario de la respuesta ------------------------------------
  Ny = Nx + Nh - 1;
  if N < Ny
    error('N = %d no alcanza: la respuesta tiene %d elementos no nulos', N, Ny);
  end

  % --- Vectores de igual dimensión N, completando con la zona nula --------
  xN = [x, zeros(1, N - Nx)];
  hN = [h, zeros(1, N - Nh)];

  % --- Convolución por FFT ------------------------------------------------
  X = fft(xN);
  H = fft(hN);
  Y = X .* H;
  yN = real(ifft(Y));          % la parte imaginaria es error de redondeo

  % El primer elemento corresponde a n = nx(1) + nh(1) (igual que conv)
  n0    = nx(1) + nh(1);
  nN    = n0 + (0:N-1);        % eje [n] de los N elementos
  ny    = n0 + (0:Ny-1);
  y_fft = yN(1:Ny);

  % --- Comparación con la convolución directa del Ejercicio 2 -------------
  y_conv = conv(x, h);
  printf('Max |y_fft - y_conv| = %.3e\n', max(abs(y_fft - y_conv)));

  % --- Gráficos: 3 filas x 2 columnas -------------------------------------
  k   = 0:N-1;                 % frecuencias discretas de la FFT
  lim = [nN(1), nN(end)];

  figure('Name', sprintf('Ejercicio 3 - Convolucion por FFT (N = %d)', N));

  subplot(3, 2, 1);
  stem(nx(1) + (0:N-1), xN, 'filled', 'markersize', 2);
  xlim(lim); grid on;
  title('x[n] completada con ceros'); ylabel('x[n]');

  subplot(3, 2, 2);
  stem(k, abs(X), 'filled', 'markersize', 2);
  xlim([0 N-1]); grid on;
  title('|X[k]|');

  subplot(3, 2, 3);
  stem(nh(1) + (0:N-1), hN, 'filled', 'markersize', 2);
  xlim(lim); grid on;
  title('h[n] completada con ceros'); ylabel('h[n]');

  subplot(3, 2, 4);
  stem(k, abs(H), 'filled', 'markersize', 2);
  xlim([0 N-1]); grid on;
  title('|H[k]|  (fenómeno de Gibbs en el corte)');

  subplot(3, 2, 5);
  stem(nN, yN, 'filled', 'markersize', 2);
  xlim(lim); grid on;
  title('y[n] = IFFT(X[k] H[k])'); xlabel('n'); ylabel('y[n]');

  subplot(3, 2, 6);
  stem(k, abs(Y), 'filled', 'markersize', 2);
  xlim([0 N-1]); grid on;
  title('|Y[k]| = |X[k]| |H[k]|'); xlabel('k');

end
