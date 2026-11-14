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
%
% QUÉ PIDE LA CONSIGNA
%   - Resolver el Ej. 2 pero haciendo la convolución con la FFT.
%   - Llevar x[n] y h[n] a la misma longitud N = 2^r (ej. 1024) agregando
%     ceros, con N suficiente para que entren todos los valores de y[n].
%   - Figura de 3 filas x 2 columnas: izquierda las señales en [n],
%     derecha el MÓDULO (abs) de sus FFT en función de [k].
%   - Comparar con el Ej. 2, notar el fenómeno de Gibbs, comparar los
%     espectros de entrada y salida.
%
% TEORÍA
%   * Propiedad de convolución: convolucionar en el tiempo equivale a
%     MULTIPLICAR en frecuencia. En vez de hacer la suma de la convolución
%     (lenta: ~Nx·Nh productos), se hace:
%         X = FFT(x),  H = FFT(h),  Y = X·H (elemento a elemento),
%         y = IFFT(Y)  (transformada inversa)
%   * La FFT es un algoritmo rápido para calcular la DFT (transformada
%     discreta de Fourier). Es más rápida cuando N es potencia de 2, por
%     eso la consigna pide N = 2^r.
%   * ¿POR QUÉ HACEN FALTA LOS CEROS? La DFT trata las señales como si
%     fueran PERIÓDICAS de periodo N. El producto X·H corresponde a una
%     convolución "circular": lo que se pasa del final vuelve a entrar por
%     el principio. Si N >= Nx + Nh - 1 (el largo de y), no hay nada que
%     se pase, y la convolución circular coincide con la lineal (la de conv).
%     Acá: Nx + Nh - 1 = 401 + 101 - 1 = 501 <= 1024. OK.
%   * QUÉ SIGNIFICA k: el índice k = 0, 1, ..., N-1 de la FFT es una
%     frecuencia discreta. Corresponde a la frecuencia analógica
%         f_k = k / (N·Δt) = k · Fm / N
%     Con Fm = 10000 Hz y N = 1024, cada k equivale a ~9.77 Hz, entonces:
%         150 Hz -> k ≈ 15,   500 Hz (Fc) -> k ≈ 51,   1000 Hz -> k ≈ 102
%     y los valores de k mayores a N/2 = 512 son la "imagen espejo" de las
%     frecuencias negativas (por eso el espectro se ve simétrico).
% =========================================================================
function [ny, y_fft, Y] = ej3(dt, Nh, Fc, Nx, N)

  % Valores por defecto (misma lógica que en ej2.m).
  if nargin < 1 || isempty(dt), dt = 0.0001; end
  if nargin < 2 || isempty(Nh), Nh = 101;    end
  if nargin < 3 || isempty(Fc), Fc = 500;    end
  if nargin < 4 || isempty(Nx), Nx = 401;    end
  if nargin < 5 || isempty(N),  N  = 1024;   end   % N = 2^r

  % --- Señales del Ejercicio 2 --------------------------------------------
  % Exactamente las mismas x[n] y h[n] del Ej. 2 (acá t = nx*dt se escribe
  % directamente dentro del cos y el sin).
  nx = 0:Nx-1;
  x  = 10*cos(2*pi*150*nx*dt) + 5*sin(2*pi*1000*nx*dt);

  M  = (Nh - 1) / 2;
  nh = -M:M;
  h  = 2*Fc*dt * sinc(2*Fc*dt*nh);

  % --- Largo necesario de la respuesta ------------------------------------
  % Ny = cantidad de valores no nulos de y (lo mismo que devolvería conv).
  % Si N no alcanza, se corta con un mensaje de error en vez de dar un
  % resultado incorrecto (convolución circular "solapada").
  Ny = Nx + Nh - 1;
  if N < Ny
    error('N = %d no alcanza: la respuesta tiene %d elementos no nulos', N, Ny);
  end

  % --- Vectores de igual dimensión N, completando con la zona nula --------
  % zeros(1, k) crea un vector fila de k ceros.
  % [x, zeros(...)] pega los ceros al final de x -> vector de largo N.
  % xN = [x(1) ... x(401) 0 0 ... 0]   (623 ceros)
  % hN = [h(1) ... h(101) 0 0 ... 0]   (923 ceros)
  % Notar que h(1) es el valor en n = -50, no en n = 0. Ese corrimiento se
  % compensa más abajo con n0.
  xN = [x, zeros(1, N - Nx)];
  hN = [h, zeros(1, N - Nh)];

  % --- Convolución por FFT ------------------------------------------------
  % fft(v) devuelve un vector de N números COMPLEJOS (módulo y fase).
  X = fft(xN);
  H = fft(hN);
  % ".*" es la multiplicación ELEMENTO A ELEMENTO: Y(k) = X(k)·H(k).
  % (Sin el punto, "*" sería producto matricial y daría error o algo
  %  distinto. Lo mismo pasa con "./" y ".^".)
  Y = X .* H;
  % ifft = FFT inversa: vuelve al dominio del tiempo. Por errores de
  % redondeo de la computadora queda una parte imaginaria minúscula
  % (~1e-15), que se descarta con real().
  yN = real(ifft(Y));          % la parte imaginaria es error de redondeo

  % El primer elemento corresponde a n = nx(1) + nh(1) (igual que conv)
  % n0 = 0 + (-50) = -50. Los N valores de yN van de n = -50 a n = 973,
  % pero solo los primeros Ny = 501 (n = -50 ... 450) son la respuesta;
  % el resto son los ceros que agregamos.
  n0    = nx(1) + nh(1);
  nN    = n0 + (0:N-1);        % eje [n] de los N elementos
  ny    = n0 + (0:Ny-1);
  % v(1:Ny) toma los elementos del 1 al Ny.
  y_fft = yN(1:Ny);

  % --- Comparación con la convolución directa del Ejercicio 2 -------------
  % Se calcula también con conv y se imprime la máxima diferencia.
  % RESULTADO: Max |y_fft - y_conv| = 1.243e-14, o sea ~0: son iguales
  % salvo errores de redondeo de la computadora (que trabaja con ~16 cifras).
  y_conv = conv(x, h);
  printf('Max |y_fft - y_conv| = %.3e\n', max(abs(y_fft - y_conv)));

  % --- Gráficos: 3 filas x 2 columnas -------------------------------------
  % subplot(3, 2, p): la posición p se cuenta por filas, de izquierda a
  % derecha:   1 2
  %            3 4
  %            5 6
  % Columna izquierda (1, 3, 5) = tiempo; derecha (2, 4, 6) = frecuencia.
  k   = 0:N-1;                 % frecuencias discretas de la FFT
  lim = [nN(1), nN(end)];

  figure('Name', sprintf('Ejercicio 3 - Convolucion por FFT (N = %d)', N));

  subplot(3, 2, 1);
  % x completada con ceros, dibujada en su verdadero eje n (0 ... 1023).
  stem(nx(1) + (0:N-1), xN, 'filled', 'markersize', 2);
  xlim(lim); grid on;
  title('x[n] completada con ceros'); ylabel('x[n]');

  subplot(3, 2, 2);
  % abs() = módulo del número complejo. Se ven DOS picos (y sus espejos
  % del otro lado de k = 512): uno alto en k ≈ 15 (150 Hz, amplitud 10)
  % y otro más bajo en k ≈ 102 (1000 Hz, amplitud 5). La altura teórica
  % del pico es amplitud·Nx/2 (2005 y 1002); da algo menos (≈ 1922 y
  % ≈ 962) porque 150 y 1000 Hz no caen exactamente sobre un k entero.
  % Los picos no son rayas finas sino que tienen "faldas" (lóbulos
  % laterales): es el efecto de haber cortado x(t) a una ventana finita
  % (enventanado rectangular, "leakage" o fuga espectral).
  stem(k, abs(X), 'filled', 'markersize', 2);
  xlim([0 N-1]); grid on;
  title('|X[k]|');

  subplot(3, 2, 3);
  % h completada con ceros, desde n = -50.
  stem(nh(1) + (0:N-1), hN, 'filled', 'markersize', 2);
  xlim(lim); grid on;
  title('h[n] completada con ceros'); ylabel('h[n]');

  subplot(3, 2, 4);
  % |H[k]| debería ser un "escalón": 1 para k < 51 (f < 500 Hz) y 0 por
  % encima (más su espejo). Pero cerca del corte aparecen OSCILACIONES y
  % un sobrepico de ~10 % por encima de 1: es el FENÓMENO DE GIBBS,
  % causado por cortar la sinc infinita a 101 valores (ventana rectangular).
  % Agrandar Nh NO lo elimina: medido con Nh = 101, 201 y 401 da 10.0 %,
  % 9.5 % y 9.2 % (tiende a ~9 %); las oscilaciones solo se angostan y se
  % pegan al corte. Para reducirlo hace falta otra ventana -> Ej. 4 (Hamming).
  stem(k, abs(H), 'filled', 'markersize', 2);
  xlim([0 N-1]); grid on;
  title('|H[k]|  (fenómeno de Gibbs en el corte)');

  subplot(3, 2, 5);
  % Resultado de la IFFT: los primeros 501 valores son idénticos a los del
  % Ej. 2, y después solo hay ceros (prueba de que N alcanzó).
  stem(nN, yN, 'filled', 'markersize', 2);
  xlim(lim); grid on;
  title('y[n] = IFFT(X[k] H[k])'); xlabel('n'); ylabel('y[n]');

  subplot(3, 2, 6);
  % Espectro de la salida: |Y| = |X|·|H|. Comparado con |X[k]|, el pico
  % de 150 Hz quedó igual (|H| ≈ 1 ahí) y el de 1000 Hz desapareció
  % (|H| ≈ 0 ahí). Esto es "ver" el filtrado en frecuencia.
  stem(k, abs(Y), 'filled', 'markersize', 2);
  xlim([0 N-1]); grid on;
  title('|Y[k]| = |X[k]| |H[k]|'); xlabel('k');

end
