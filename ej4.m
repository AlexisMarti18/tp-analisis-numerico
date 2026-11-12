% =========================================================================
% ANÁLISIS NUMÉRICO - TRABAJOS DE LABORATORIO
% EJERCICIO 4: Filtro pasa bajos FIR por sinc enventanada (Hamming)
%
% La frecuencia de muestreo Fm es el dato principal: dt = 1/Fm.
%   h[n] = (2 Fc/Fm) sinc(2 Fc n / Fm) . w_Hamming[n],  normalizado a sum(h) = 1
%
% Uso:  ej4                         -> valores por defecto
%       ej4(Fm, Nh, Fc, Nx, N)
% =========================================================================
function [t_y, y_conv, y_fft, f, X, H, Y] = ej4(Fm, Nh, Fc, Nx, N)

  if nargin < 1 || isempty(Fm), Fm = 10000; end  % frecuencia de muestreo [Hz]
  if nargin < 2 || isempty(Nh), Nh = 101;   end  % largo del núcleo (impar)
  if nargin < 3 || isempty(Fc), Fc = 500;   end  % frecuencia de corte [Hz]
  if nargin < 4 || isempty(Nx), Nx = 401;   end  % muestras de la entrada
  if nargin < 5 || isempty(N),  N  = max(1024, 2^nextpow2(Nx + Nh - 1)); end

  if mod(Nh, 2) == 0, error('Nh debe ser impar'); end
  if N < Nx + Nh - 1, error('N debe ser >= Nx + Nh - 1'); end

  % --- Entrada muestreada a Fm ----------------------------------------------
  nx  = 0:Nx-1;
  t_x = nx / Fm;                                   % tiempos continuos [s]
  x   = 10*cos(2*pi*150*t_x) + 5*sin(2*pi*1000*t_x);

  % --- Núcleo: sinc truncada x ventana de Hamming, normalizada --------------
  M   = (Nh - 1) / 2;
  nh  = -M:M;
  t_h = nh / Fm;
  h   = (2*Fc/Fm) * sinc(2*Fc*nh/Fm);              % sinc ideal truncada
  h   = h .* hamming(Nh)';                         % ventana de Hamming
  h   = h / sum(h);                                % normalizado: H(0) = 1

  % --- Respuesta por convolución en el tiempo -------------------------------
  y_conv = conv(x, h);
  ny     = (nx(1) + nh(1)) : (nx(end) + nh(end));
  t_y    = ny / Fm;

  % --- Respuesta por convolución FFT ----------------------------------------
  X     = fft(x, N);                               % fft(x, N) completa con ceros
  H     = fft(h, N);
  Y     = X .* H;
  yN    = real(ifft(Y));
  y_fft = yN(1:numel(ny));

  printf('Fm = %g Hz | Nh = %d | N = %d | max|y_fft - y_conv| = %.3e\n', ...
         Fm, Nh, N, max(abs(y_fft - y_conv)));

  % --- Eje de frecuencias analógicas de la FFT ------------------------------
  f    = (0:N-1) * Fm / N;                         % f_k = k Fm / N  [Hz]
  flim = [0, Fm/2];                                % hasta Nyquist (espectro simétrico)
  tlim = [t_y(1), t_y(end)];

  % --- Gráficos: 4 filas x 2 columnas ---------------------------------------
  figure('Name', sprintf('Ejercicio 4 - FIR Hamming (Fm = %g Hz, Nh = %d)', Fm, Nh));

  subplot(4, 2, 1);
  plot(t_x, x, '.-'); xlim(tlim); grid on;
  title(sprintf('Entrada x(t_n)   (F_m = %g Hz, N_x = %d)', Fm, Nx));
  ylabel('x');

  subplot(4, 2, 2);
  plot(f, abs(X)); xlim(flim); grid on;
  title('|X(f)|');

  subplot(4, 2, 3);
  plot(t_h, h, '.-'); xlim(tlim); grid on;
  title(sprintf('h(t_n): sinc x Hamming normalizada  (N_h = %d, F_c = %g Hz)', Nh, Fc));
  ylabel('h');

  subplot(4, 2, 4);
  plot(f, abs(H)); xlim(flim); grid on;
  title('|H(f)|');

  subplot(4, 2, 5);
  plot(t_y, y_conv, '.-'); xlim(tlim); grid on;
  title('Salida y(t_n) por convolución en el tiempo (conv)');
  ylabel('y');

  subplot(4, 2, 6);
  plot(f, abs(Y)); xlim(flim); grid on;
  title('|Y(f)| = |X(f)| |H(f)|  (FFT)');

  subplot(4, 2, 7);
  plot(t_y, y_fft, '.-'); xlim(tlim); grid on;
  title('Salida y(t_n) por convolución FFT (ifft)');
  xlabel('t [s]'); ylabel('y');

  subplot(4, 2, 8);
  plot(f, abs(Y)); xlim(flim); grid on;
  title('|Y(f)| = |X(f)| |H(f)|  (FFT)');
  xlabel('f [Hz]');

end
