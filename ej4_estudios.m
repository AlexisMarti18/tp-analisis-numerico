% =========================================================================
% EJERCICIO 4 (continuación): estudios pedidos en el enunciado
%   a) Mejora de H(f) con Hamming frente a la ventana rectangular (Gibbs)
%   b) Banda de transición al modificar el núcleo Nh
%   c) Espectro de la entrada al modificar la ventana no nula Nx
%   d) Calidad de representación al modificar la frecuencia de muestreo Fm
% =========================================================================
%
% El final de la consigna del Ej. 4 pide "verificar" e "investigar" cuatro
% cosas. Este script hace un estudio por cada una, con su propia figura y
% números impresos en consola para poder fundamentar las conclusiones.
%
% SINTAXIS NUEVA EN ESTE ARCHIVO
%   dB = @(v) 20*log10(abs(v) + 1e-12);
%       Función ANÓNIMA: define una mini-función en una línea, sin archivo
%       aparte. Después se usa como dB(H). Pasa un módulo a decibeles:
%       1 -> 0 dB, 0.1 -> -20 dB, 0.01 -> -40 dB, 0.001 -> -60 dB.
%       El + 1e-12 evita calcular log10(0) = -infinito.
%   plot(x1, y1, x2, y2)  -> dos curvas en el mismo gráfico.
%   hold on / hold off    -> mientras está "on", cada plot nuevo se agrega
%                            al gráfico en vez de borrar el anterior.
%   legend('a', 'b')      -> recuadro con el nombre de cada curva.
%   v(condición)          -> "indexado lógico": se queda solo con los
%                            elementos de v donde la condición es verdadera.
%                            Ej.: H(f > 650 & f < 5000)  ("&" = Y lógico).
%   find(cond, 1)         -> posición del PRIMER elemento que cumple cond.
%   numel(v)              -> cantidad de elementos de v.
%   ~                     -> en [a, b, ~, d] = f(...), descarta esa salida.
% =========================================================================
clear; clc; close all;

% Parámetros comunes. Nfft = 8192 es una FFT mucho más larga que la señal:
% completar con muchos ceros NO agrega información, pero calcula el
% espectro en muchos más puntos (cada 10000/8192 ≈ 1.2 Hz en vez de
% ≈ 9.8 Hz), así la curva se ve continua y se pueden medir bien los
% anchos y los lóbulos.
Fc = 500;  Fm = 10000;  Nfft = 8192;      % FFT larga: H(f) y X(f) casi continuos
f  = (0:Nfft-1) * Fm / Nfft;
dB = @(v) 20*log10(abs(v) + 1e-12);

% ---------- a) Rectangular vs Hamming -------------------------------------
% Consigna: "Verificar la mejora en el espectro H(f) reduciendo el
% fenómeno de Gibbs". Se calculan los dos núcleos con el mismo Nh = 101:
%   h_rect = sinc truncada (la de los Ej. 2 y 3)
%   h_hamm = sinc x Hamming normalizada (la del Ej. 4)
Nh = 101;  nh = -(Nh-1)/2:(Nh-1)/2;
h_rect = (2*Fc/Fm) * sinc(2*Fc*nh/Fm);
h_hamm = h_rect .* hamming(Nh)';   h_hamm = h_hamm / sum(h_hamm);
H_rect = fft(h_rect, Nfft);        H_hamm = fft(h_hamm, Nfft);

% Figura con 2 gráficos: |H| en escala lineal y en decibeles. La escala en
% dB permite ver los lóbulos laterales chiquitos (en lineal parecen 0).
figure('Name', 'a) Rectangular vs Hamming');
subplot(2, 1, 1);
plot(f, abs(H_rect), f, abs(H_hamm)); xlim([0 Fm/2]); grid on;
legend('Rectangular (Ej. 2 y 3)', 'Hamming normalizada (Ej. 4)');
title(sprintf('|H(f)|, N_h = %d, F_c = %d Hz', Nh, Fc));
subplot(2, 1, 2);
plot(f, dB(H_rect), f, dB(H_hamm)); xlim([0 Fm/2]); ylim([-100 5]); grid on;
title('|H(f)| en dB'); xlabel('f [Hz]'); ylabel('dB');

% Mediciones:
%   sobrepico = cuánto supera max|H| al valor ideal 1, en %.
%   lóbulo lateral máximo = el mayor |H| en la banda eliminada, tomada
%   como f entre 1.3·Fc (650 Hz) y Nyquist, en dB.
% RESULTADOS:
%   Rectangular: sobrepico = 10.0 %, lóbulo lateral máx = -27.9 dB
%   Hamming    : sobrepico =  0.0 %, lóbulo lateral máx = -41.0 dB
% -> Hamming elimina el sobrepico y deja pasar ~4.5 veces menos amplitud
%    en la banda eliminada (13 dB menos). Gibbs reducido: verificado.
%    (A 650 Hz Hamming todavía está terminando su transición; más lejos
%     del corte sus lóbulos son aún más bajos, se ve en la figura en dB.)
printf('Rectangular: sobrepico = %.1f %%, lobulo lateral max = %.1f dB\n', ...
       (max(abs(H_rect)) - 1)*100, max(dB(H_rect(f > 1.3*Fc & f < Fm/2))));
printf('Hamming    : sobrepico = %.1f %%, lobulo lateral max = %.1f dB\n', ...
       (max(abs(H_hamm)) - 1)*100, max(dB(H_hamm(f > 1.3*Fc & f < Fm/2))));

% ---------- b) Banda de transición vs Nh -----------------------------------
% Consigna: "Investigar los cambios en la calidad de su comportamiento en
% la banda de transición, donde el filtro ideal es discontinuo, al
% modificar las dimensiones del núcleo del filtro (h[n])".
% Se dibuja |H(f)| (con Hamming) para 4 valores de Nh en el mismo gráfico.
figure('Name', 'b) Banda de transicion vs Nh');
hold on;
lista_Nh = [21 51 101 201];
for Nh = lista_Nh
  nh = -(Nh-1)/2:(Nh-1)/2;
  h  = (2*Fc/Fm) * sinc(2*Fc*nh/Fm) .* hamming(Nh)';  h = h / sum(h);
  H  = abs(fft(h, Nfft));
  plot(f, H);
  % ancho de transición: de |H| = 0.9 a |H| = 0.1
  % f1 = primera frecuencia donde |H| baja de 0.9 (fin de la banda de paso)
  % f2 = primera frecuencia donde |H| baja de 0.1 (inicio de la eliminada)
  % También se mide cuánto deja pasar en 1000 Hz (la componente a eliminar).
  % round(1000*Nfft/Fm) + 1 = índice del vector f más cercano a 1000 Hz
  % (el +1 es porque en Octave los índices empiezan en 1 y f(1) = 0 Hz).
  f1 = f(find(H < 0.9, 1));  f2 = f(find(H < 0.1, 1));
  printf('Nh = %3d: transicion %4.0f - %4.0f Hz (ancho %4.0f Hz), |H(1000 Hz)| = %.4f\n', ...
         Nh, f1, f2, f2 - f1, H(round(1000*Nfft/Fm) + 1));
end
% RESULTADOS:
%   Nh =  21: transición 228 - 989 Hz (ancho 760 Hz), |H(1000 Hz)| = 0.0933
%   Nh =  51: transición 316 - 681 Hz (ancho 365 Hz), |H(1000 Hz)| = 0.0017
%   Nh = 101: transición 409 - 592 Hz (ancho 183 Hz), |H(1000 Hz)| = 0.0007
%   Nh = 201: transición 455 - 546 Hz (ancho  90 Hz), |H(1000 Hz)| = 0.0005
% -> El ancho de la transición es INVERSAMENTE PROPORCIONAL a Nh: duplicar
%    Nh lo reduce a la mitad (Nh · ancho ≈ 18000 ≈ 1.8·Fm, constante).
%    Con Nh = 21 el filtro es tan "blando" que en 1000 Hz todavía deja
%    pasar 9 % (0.47 de amplitud de los 5 originales). Con Nh >= 51 ya es
%    despreciable. Más Nh = filtro más parecido al ideal, pero más cuentas
%    y más retardo/transitorio.
hold off; xlim([0 2000]); grid on;
% arrayfun aplica una función a cada elemento de lista_Nh; acá arma los
% textos 'N_h = 21', 'N_h = 51', ... para la leyenda.
legend(arrayfun(@(v) sprintf('N_h = %d', v), lista_Nh, 'UniformOutput', false));
title('|H(f)| con ventana de Hamming para distintos N_h'); xlabel('f [Hz]');

% ---------- c) Espectro de la entrada vs ventana no nula Nx ---------------
% Consigna: "Investigar los cambios en la calidad del espectro de la señal
% de entrada al modificar el tamaño de la ventana no nula, y verificar el
% efecto enventanado comparándolo con el espectro de la señal de longitud
% infinita".
%   - Una señal de longitud INFINITA con dos senoidales tiene un espectro
%     de dos rayas finas (deltas): altura 10 en 150 Hz y 5 en 1000 Hz.
%     Se dibujan en rojo con stem como referencia.
%   - Tomar solo Nx muestras = multiplicar por una ventana rectangular.
%     En frecuencia, cada raya se convierte en un "lóbulo" de ancho
%     ~2·Fm/Nx más lóbulos laterales (fuga espectral o leakage).
figure('Name', 'c) Efecto de enventanado sobre X(f)');
lista_Nx = [101 401 1601];
% "for i = 1:numel(lista)" recorre las posiciones 1, 2, 3 para poder
% usar i también como número de subplot.
for i = 1:numel(lista_Nx)
  Nx = lista_Nx(i);  n = 0:Nx-1;
  x  = 10*cos(2*pi*150*n/Fm) + 5*sin(2*pi*1000*n/Fm);
  % Escalado 2|X|/Nx: sin esto la altura del pico crece con Nx. Así, el
  % pico de una senoidal de amplitud A tiende a A, comparable con las rayas.
  X  = abs(fft(x, Nfft)) / Nx * 2;       % escalado: el pico tiende a la amplitud
  subplot(numel(lista_Nx), 1, i);
  plot(f, X); hold on;
  stem([150 1000], [10 5], 'r', 'filled'); hold off;   % espectro de x(t) infinita
  xlim([0 1500]); ylim([0 11]); grid on;
  title(sprintf('2|X(f)|/N_x con N_x = %d (ventana de %.1f ms)', Nx, 1000*Nx/Fm));
end
% QUÉ SE OBSERVA
%   Nx =  101 (10 ms):  lóbulos anchos (~200 Hz), los picos se "desparraman".
%   Nx =  401 (40 ms):  lóbulos de ~50 Hz.
%   Nx = 1601 (160 ms): lóbulos de ~12 Hz, muy parecido a las dos rayas.
% -> Cuanto más larga la ventana, más se parece al espectro de la señal
%    infinita (mejor RESOLUCIÓN en frecuencia). Con una ventana muy corta,
%    dos frecuencias cercanas se mezclarían en un solo lóbulo.
xlabel('f [Hz]');

% ---------- d) Variación de la frecuencia de muestreo Fm ------------------
% Consigna: "Verificar las modificaciones en la calidad de representación
% de las señales en tiempo discreto y frecuencias, al modificar Fm".
% Se usa la función ej4 con 4 valores de Fm, manteniendo FIJAS en segundos
% la duración de la ventana (T = 40 ms) y la del núcleo (~10 ms), para que
% lo único que cambie sea la densidad de muestras.
fd = figure('Name', 'd) Variacion de Fm');   % fd guarda el "identificador" de la figura
lista_Fm = [10000 3000 2000 1200];
T = 0.04;                                 % misma duración de ventana [s]
for i = 1:numel(lista_Fm)
  Fmi = lista_Fm(i);
  Nx  = round(T*Fmi) + 1;                 % muestras en 40 ms: 401, 121, 81, 49
  Nh  = 2*round(0.005*Fmi) + 1;           % núcleo de ~10 ms en todos los casos
  % Se llama a ej4, que abre su propia figura de 4x2; se usan sus resultados
  % y se cierra esa figura (close(gcf); gcf = "figura actual").
  [t_y, y_conv, ~, fi, X] = ej4(Fmi, Nh, Fc, Nx);  close(gcf);
  figure(fd);                             % volver a la figura del estudio
  t_x = (0:Nx-1) / Fmi;
  x   = 10*cos(2*pi*150*t_x) + 5*sin(2*pi*1000*t_x);

  subplot(numel(lista_Fm), 2, 2*i-1);     % columna izquierda: tiempo
  plot(t_x, x, '.-', t_y, y_conv, '.-'); xlim([0 T]); grid on;
  title(sprintf('F_m = %d Hz: x (azul) e y (rojo)', Fmi));
  subplot(numel(lista_Fm), 2, 2*i);       % columna derecha: espectro
  plot(fi, abs(X)); xlim([0 Fmi/2]); grid on;
  title(sprintf('|X(f)|, F_m = %d Hz (Nyquist = %d Hz)', Fmi, Fmi/2));
end
% QUÉ SE OBSERVA (regla: se representa bien lo que está por debajo de Fm/2)
%   Fm = 10000 Hz (Nyquist 5000): las dos componentes bien representadas;
%     x se ve suave en el tiempo y el espectro tiene sus dos picos.
%   Fm =  3000 Hz (Nyquist 1500): 1000 Hz < 1500, se representa, pero con
%     solo 3 muestras por periodo: en el tiempo x se ve "quebrada". El
%     espectro igual muestra los dos picos en su lugar y el filtro funciona.
%   Fm =  2000 Hz (Nyquist 1000): la componente de 1000 Hz cae justo en
%     Nyquist y, al ser un seno, todas sus muestras valen 0: desaparece.
%   Fm =  1200 Hz (Nyquist 600): ALIASING. La componente de 1000 Hz aparece
%     disfrazada en 1200 - 1000 = 200 Hz. Como 200 < Fc = 500 Hz, el filtro
%     la deja PASAR y la salida queda contaminada. Ningún filtro digital
%     puede arreglar esto después de muestrear: hay que cumplir Nyquist
%     (o filtrar en analógico ANTES de muestrear, filtro "antialiasing").
xlabel('f [Hz]');
