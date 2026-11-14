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
%
% QUÉ PIDE LA CONSIGNA
%   - Implementar un filtro pasa bajos con la misma señal y sistema.
%   - Usar la frecuencia de muestreo Fm como dato, en lugar de Δt (Δt = 1/Fm).
%   - Multiplicar h[n] por una ventana de Hamming (comando hamming; el
%     enunciado dice "haming", que es un error de tipeo) para reducir el
%     fenómeno de Gibbs, y después normalizar.
%   - Calcular la salida por conv Y por FFT, y comparar.
%   - Figura de 4 filas x 2 columnas: izquierda señales vs TIEMPO EN
%     SEGUNDOS, derecha módulos de los espectros vs FRECUENCIA EN Hz.
%     Para el espectro de las dos salidas se usa el mismo Y de la FFT.
%   - Estudios adicionales (Nh, Nx, Fm, Gibbs) -> ver ej4_estudios.m.
%
% TEORÍA: ¿QUÉ ES UNA VENTANA Y POR QUÉ HAMMING?
%   * En los Ej. 2 y 3 la sinc se cortó "de golpe" a Nh valores. Eso es
%     multiplicarla por una ventana RECTANGULAR (1 adentro, 0 afuera). El
%     corte brusco en el tiempo produce oscilaciones en frecuencia (Gibbs):
%     sobrepico de ~9 % en la banda de paso y lóbulos laterales altos en
%     la banda eliminada (deja "escapar" frecuencias que debería anular).
%   * La ventana de Hamming es una campana suave:
%         w[m] = 0.54 - 0.46·cos(2π·m/(Nh-1)),   m = 0 ... Nh-1
%     vale 1 en el centro y ~0.08 en los extremos. Al multiplicar la sinc
%     por w, los extremos se "apagan" suavemente en vez de cortarse de golpe
%     -> desaparece el sobrepico y los lóbulos laterales bajan mucho.
%   * El precio a pagar: la banda de transición (el paso de 1 a 0 en |H|)
%     se hace más ANCHA que con la rectangular del mismo Nh.
%   * Normalizar: al aplicar la ventana, la suma de h cambia, y la suma de
%     h es justamente H(0) (la ganancia en continua / bajas frecuencias).
%     Dividiendo h por sum(h) se fuerza H(0) = 1: las frecuencias bajas
%     (como 150 Hz) pasan con su amplitud original.
%   * "FIR" = Finite Impulse Response: h[n] tiene una cantidad finita de
%     valores (Nh). Este método se llama "diseño por ventanas".
%
% POR QUÉ USAR Fm EN VEZ DE Δt
%   Fm es el parámetro con significado físico: dice qué frecuencias se
%   pueden representar (hasta Fm/2, Nyquist). Todas las fórmulas cambian
%   Δt por 1/Fm:  t = n/Fm,  h[n] = (2Fc/Fm)·sinc(2Fc·n/Fm),  f_k = k·Fm/N.
% =========================================================================
function [t_y, y_conv, y_fft, f, X, H, Y] = ej4(Fm, Nh, Fc, Nx, N)

  % --- Parámetros (Fm reemplaza a dt) --------------------------------------
  if nargin < 1 || isempty(Fm), Fm = 10000; end  % frecuencia de muestreo [Hz]
  if nargin < 2 || isempty(Nh), Nh = 101;   end  % largo del núcleo (impar)
  if nargin < 3 || isempty(Fc), Fc = 500;   end  % frecuencia de corte [Hz]
  if nargin < 4 || isempty(Nx), Nx = 401;   end  % muestras de la entrada
  % N por defecto: la menor potencia de 2 que alcance para Nx + Nh - 1,
  % pero como mínimo 1024. nextpow2(v) da el exponente r tal que 2^r >= v.
  % Ej.: Nx + Nh - 1 = 501 -> nextpow2 = 9 -> 2^9 = 512 -> max(1024, 512) = 1024.
  if nargin < 5 || isempty(N),  N  = max(1024, 2^nextpow2(Nx + Nh - 1)); end

  % Validaciones: Nh impar (h centrada) y N suficiente (convolución lineal).
  if mod(Nh, 2) == 0, error('Nh debe ser impar'); end
  if N < Nx + Nh - 1, error('N debe ser >= Nx + Nh - 1'); end

  % --- Entrada muestreada a Fm ----------------------------------------------
  % Igual que antes pero con t = n / Fm (en vez de n·Δt).
  nx  = 0:Nx-1;
  t_x = nx / Fm;                                   % tiempos continuos [s]
  x   = 10*cos(2*pi*150*t_x) + 5*sin(2*pi*1000*t_x);

  % --- Núcleo: sinc truncada x ventana de Hamming, normalizada --------------
  % Se construye en 3 pasos (cada línea REEMPLAZA a h por la nueva versión):
  M   = (Nh - 1) / 2;
  nh  = -M:M;
  t_h = nh / Fm;                                   % h centrada en t = 0
  % 1) La misma sinc truncada de los Ej. 2 y 3 (2FcΔt = 2Fc/Fm).
  h   = (2*Fc/Fm) * sinc(2*Fc*nh/Fm);              % sinc ideal truncada
  % 2) hamming(Nh) devuelve la ventana como vector COLUMNA (Nh x 1); el
  %    apóstrofo ' la transpone a FILA (1 x Nh) para que tenga la misma
  %    forma que h. Luego ".*" multiplica valor a valor: h[n]·w[n].
  h   = h .* hamming(Nh)';                         % ventana de Hamming
  % 3) Normalizado: sum(h) suma todos los valores. Después de esto
  %    sum(h) = 1, es decir H(f = 0) = 1 (ganancia unitaria en continua).
  h   = h / sum(h);                                % normalizado: H(0) = 1

  % --- Respuesta por convolución en el tiempo -------------------------------
  % Igual que en el Ej. 2; el eje n se pasa a segundos dividiendo por Fm.
  y_conv = conv(x, h);
  ny     = (nx(1) + nh(1)) : (nx(end) + nh(end));
  t_y    = ny / Fm;

  % --- Respuesta por convolución FFT ----------------------------------------
  % fft(v, N) es un atajo: completa v con ceros hasta largo N y calcula la
  % FFT. Hace lo mismo que [v, zeros(1, N - length(v))] del Ej. 3.
  X     = fft(x, N);                               % fft(x, N) completa con ceros
  H     = fft(h, N);
  Y     = X .* H;
  yN    = real(ifft(Y));
  % numel(v) = cantidad de elementos. Se toman los primeros Nx+Nh-1.
  y_fft = yN(1:numel(ny));

  % Comparación numérica de los dos métodos.
  % RESULTADO: max|y_fft - y_conv| = 1.066e-14 -> idénticos (redondeo).
  printf('Fm = %g Hz | Nh = %d | N = %d | max|y_fft - y_conv| = %.3e\n', ...
         Fm, Nh, N, max(abs(y_fft - y_conv)));

  % --- Eje de frecuencias analógicas de la FFT ------------------------------
  % El elemento k de la FFT corresponde a la frecuencia f_k = k·Fm/N [Hz].
  % Con Fm = 10000 y N = 1024: f = [0, 9.77, 19.5, ..., 9990] Hz.
  f    = (0:N-1) * Fm / N;                         % f_k = k Fm / N  [Hz]
  % Solo se grafica de 0 a Fm/2 (frecuencia de Nyquist): la mitad superior
  % del espectro es el espejo de la inferior (no aporta información nueva).
  flim = [0, Fm/2];                                % hasta Nyquist (espectro simétrico)
  tlim = [t_y(1), t_y(end)];

  % --- Gráficos: 4 filas x 2 columnas ---------------------------------------
  % Posiciones:  1 2   <- entrada x     | |X(f)|
  %              3 4   <- núcleo h      | |H(f)|
  %              5 6   <- salida (conv) | |Y(f)|
  %              7 8   <- salida (FFT)  | |Y(f)| (el mismo, como pide la consigna)
  % Acá se usa plot (curva que une los puntos) en vez de stem, porque con
  % el eje en segundos y Hz se quiere ver la señal "analógica" representada.
  % '.-' = marca un punto en cada muestra y los une con una línea.
  figure('Name', sprintf('Ejercicio 4 - FIR Hamming (Fm = %g Hz, Nh = %d)', Fm, Nh));

  subplot(4, 2, 1);
  plot(t_x, x, '.-'); xlim(tlim); grid on;
  title(sprintf('Entrada x(t_n)   (F_m = %g Hz, N_x = %d)', Fm, Nx));
  ylabel('x');

  subplot(4, 2, 2);
  % Dos picos: en 150 Hz y en 1000 Hz.
  plot(f, abs(X)); xlim(flim); grid on;
  title('|X(f)|');

  subplot(4, 2, 3);
  % h se ve como una sinc que se apaga suavemente en los bordes (por la
  % ventana), en vez de cortarse de golpe como en el Ej. 2.
  plot(t_h, h, '.-'); xlim(tlim); grid on;
  title(sprintf('h(t_n): sinc x Hamming normalizada  (N_h = %d, F_c = %g Hz)', Nh, Fc));
  ylabel('h');

  subplot(4, 2, 4);
  % |H(f)|: vale 1 hasta cerca de 500 Hz, baja suavemente y queda en ~0.
  % A diferencia del Ej. 3, SIN sobrepico ni oscilaciones (Gibbs reducido).
  plot(f, abs(H)); xlim(flim); grid on;
  title('|H(f)|');

  subplot(4, 2, 5);
  % Salida: cosenoidal de 150 Hz y amplitud ~10, con transitorios en los
  % bordes. Arranca en t negativo porque h está centrada (no causal).
  plot(t_y, y_conv, '.-'); xlim(tlim); grid on;
  title('Salida y(t_n) por convolución en el tiempo (conv)');
  ylabel('y');

  subplot(4, 2, 6);
  % Solo queda el pico de 150 Hz: el de 1000 Hz fue eliminado.
  plot(f, abs(Y)); xlim(flim); grid on;
  title('|Y(f)| = |X(f)| |H(f)|  (FFT)');

  subplot(4, 2, 7);
  % Visualmente idéntica a la de arriba (fila 3).
  plot(t_y, y_fft, '.-'); xlim(tlim); grid on;
  title('Salida y(t_n) por convolución FFT (ifft)');
  xlabel('t [s]'); ylabel('y');

  subplot(4, 2, 8);
  plot(f, abs(Y)); xlim(flim); grid on;
  title('|Y(f)| = |X(f)| |H(f)|  (FFT)');
  xlabel('f [Hz]');

end
