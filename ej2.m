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
%
% QUÉ PIDE LA CONSIGNA
%   1) Muestrear x(t) = 10 cos(2π·150·t) + 5 sen(2π·1000·t)
%      con Nx muestras separadas Δt (ej.: 401 muestras, Δt = 0.0001 s),
%      empezando en t = 0.
%   2) Armar la respuesta al impulso h[n] = 2·Fc·Δt · sinc(2·Fc·Δt·n),
%      centrada en n = 0, con Nh (impar) valores: n = -(Nh-1)/2 ... (Nh-1)/2.
%   3) y[n] = conv(x, h) y graficar las 3 señales como en el Ej. 1, con el
%      eje n limitado a los elementos no nulos de y.
%   4) (en ej2_variantes.m) Repetir con Δt = 0.0005 y 0.001, y otros Nh.
%
% TEORÍA EN 4 RENGLONES
%   - Un filtro pasa bajos IDEAL deja pasar sin cambios todo lo que está
%     por debajo de Fc y elimina todo lo que está por encima.
%   - Su respuesta al impulso en tiempo continuo es h(t) = 2Fc·sinc(2Fc·t).
%   - Muestreando en t = n·Δt y multiplicando por Δt (para que la suma
%     discreta aproxime a la integral de convolución) sale exactamente
%     h[n] = 2·Fc·Δt·sinc(2·Fc·Δt·n), que es la fórmula del enunciado.
%   - La sinc ideal es infinitamente larga; acá la CORTAMOS a Nh valores
%     (= la multiplicamos por una ventana rectangular). Ese corte es lo que
%     produce el fenómeno de Gibbs que se ve en el Ej. 3.
%   Con Fc = 500 Hz: la componente de 150 Hz pasa y la de 1000 Hz se elimina.
%   Se espera ver a la salida prácticamente solo 10 cos(2π·150·t).
%
% ¿POR QUÉ ESTO ES UNA *FUNCIÓN* Y NO UN SCRIPT COMO EL EJ. 1?
%   Porque la consigna pide resolverlo varias veces con distintos Δt y Nh.
%   Una función recibe parámetros, así no hay que copiar el código.
%   El archivo se tiene que llamar igual que la función (ej2.m <-> ej2).
%
% SINTAXIS DE UNA FUNCIÓN EN OCTAVE
%   function [salida1, salida2, ...] = nombre(entrada1, entrada2, ...)
%      ... cuerpo ...
%   end
%   - Las variables de adentro son locales: no se ven desde afuera.
%   - Se puede llamar con menos entradas de las declaradas; las que faltan
%     quedan sin definir, y por eso abajo se les da un valor por defecto.
%   - Se puede pedir menos salidas: [ny, y] = ej2(...) toma solo las 2 primeras.
% =========================================================================
function [ny, y, nx, x, nh, h] = ej2(dt, Nh, Fc, Nx)

  % --- Valores por defecto de los parámetros ------------------------------
  % nargin = "number of arguments in" = cuántos parámetros se pasaron.
  % isempty(v) es verdadero si v = [] (vacío). Eso permite llamar, por
  % ejemplo, ej2([], 51) para cambiar solo Nh y dejar dt por defecto.
  % "||" es el O lógico (si la primera condición ya es verdadera, ni
  % evalúa la segunda; por eso no falla si dt no existe).
  % Forma compacta "if condición, acción; end" en una sola línea.
  if nargin < 1 || isempty(dt), dt = 0.0001; end  % paso de muestreo [s]
  if nargin < 2 || isempty(Nh), Nh = 101;    end  % largo del núcleo (impar)
  % Fc = 500 Hz (y no 1000 como el ejemplo del enunciado) porque la consigna
  % pide que Fc esté ENTRE 150 y 1000 Hz. Con 1000 Hz la componente alta
  % cae justo en el corte (ver caso c de ej2_variantes.m).
  if nargin < 3 || isempty(Fc), Fc = 500;    end  % corte [Hz], 150 < Fc < 1000
  if nargin < 4 || isempty(Nx), Nx = 401;    end  % muestras de x[n]

  % mod(a, b) = resto de dividir a por b. Si Nh es par, el resto es 0 y
  % no habría una muestra en n = 0 ni simetría -> se corta con error().
  if mod(Nh, 2) == 0
    error('Nh debe ser impar para que h[n] quede centrado en n = 0');
  end

  % --- Entrada: muestras de x(t) = 10cos(2pi150t) + 5sen(2pi1000t) -------
  % nx = [0 1 2 ... 400]: índices de muestra, empezando en t = 0.
  nx = 0:Nx-1;                 % inicio en t = 0
  % Tiempo real de cada muestra: t_n = n·Δt. Multiplicar un vector por un
  % número multiplica cada elemento: t = [0, 0.0001, 0.0002, ..., 0.04].
  t  = nx * dt;
  % cos, sin y pi existen en Octave. Al aplicarlas a un vector, se aplican
  % elemento a elemento, así que x queda con los 401 valores muestreados.
  % (Con 401 muestras de 0.0001 s, la ventana dura 0.04 s = 6 periodos de
  % la componente de 150 Hz y 40 periodos de la de 1000 Hz.)
  x  = 10*cos(2*pi*150*t) + 5*sin(2*pi*1000*t);

  % --- Respuesta al impulso: h[n] = 2 Fc dt sinc(2 Fc dt n), centrada -----
  % M = mitad del largo. Con Nh = 101: M = 50, nh = [-50 ... 50].
  M  = (Nh - 1) / 2;
  nh = -M:M;                   % de -(Nh-1)/2 a (Nh-1)/2
  % sinc en Octave está definida como sinc(u) = sin(π·u)/(π·u), con
  % sinc(0) = 1. Es la definición "normalizada", la que usa el enunciado.
  % El producto 2*Fc*dt*nh es un número por un vector -> otro vector.
  % h es simétrica (h[-n] = h[n]) y su valor máximo es h[0] = 2·Fc·Δt.
  % Con Fc = 500 y Δt = 0.0001: h[0] = 0.1, y cruza por cero cada
  % 1/(2·Fc·Δt) = 10 muestras.
  h  = 2*Fc*dt * sinc(2*Fc*dt*nh);

  % --- Respuesta del sistema ----------------------------------------------
  % Igual que en el Ej. 1: conv da los valores, y el eje n se calcula con
  % la regla "inicio = nx(1)+nh(1), fin = nx(end)+nh(end)".
  % Acá: inicio = 0 + (-50) = -50, fin = 400 + 50 = 450 -> 501 valores
  % (= Nx + Nh - 1 = 401 + 101 - 1).
  % Que y empiece en n negativo es consecuencia de que h está centrada
  % (filtro "no causal": usa muestras futuras). La ventaja es que la salida
  % NO queda retrasada respecto de la entrada.
  y  = conv(x, h);
  ny = (nx(1) + nh(1)) : (nx(end) + nh(end));

  % --- Gráficos: misma escala, límites = elementos no nulos de y[n] -------
  % La consigna pide que los límites coincidan con los elementos no nulos
  % de la respuesta, así que acá no se agrega margen (a diferencia del Ej. 1).
  lim = [ny(1), ny(end)];

  % sprintf arma un texto reemplazando los "huecos" por valores:
  %   %g -> número en el formato más corto, %d -> entero.
  % Así el nombre de la ventana dice con qué parámetros se generó.
  figure('Name', sprintf('Ejercicio 2 - dt = %g s, Nh = %d, Fc = %g Hz', dt, Nh, Fc));

  subplot(3, 1, 1);
  % 'markersize', 2 achica los círculos, porque con cientos de muestras
  % los círculos normales tapan todo.
  stem(nx, x, 'filled', 'markersize', 2);
  xlim(lim); grid on;
  % En los títulos, "_" pone subíndice y "\Delta" la letra Δ (formato TeX).
  % Dentro de sprintf, la barra "\" se escribe doble "\\" para que llegue
  % una sola al título. 1/dt es la frecuencia de muestreo Fm.
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
  % QUÉ SE VE EN y[n]:
  %  - En el centro: una cosenoidal "limpia" de amplitud ~10 y 150 Hz
  %    (la componente de 1000 Hz desapareció).
  %  - En los bordes (las primeras y últimas ~Nh muestras): un TRANSITORIO,
  %    la salida sube/baja gradualmente porque h solo se superpone en parte
  %    con x (fuera de la ventana x vale 0).
  title('Respuesta y[n] = x[n] * h[n]');
  xlabel('n'); ylabel('y[n]');

end
