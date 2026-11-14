% =========================================================================
% ANÁLISIS NUMÉRICO - TRABAJOS DE LABORATORIO
% EJERCICIO 1: Convolución discreta de un pulso rectangular con una rampa
% =========================================================================
%
% QUÉ PIDE LA CONSIGNA
%   - x[n]: pulso rectangular de valor 1, definido sobre un vector nx que
%     tenga valores de n negativos y positivos.
%   - h[n] = n, definido sobre un vector nh que empiece en n = 0.
%   - Calcular y[n] = x[n] * h[n] (convolución) con el comando conv.
%   - Graficar las 3 señales una debajo de la otra (subplot + stem), todas
%     con el mismo eje horizontal (xlim) que abarque toda la convolución.
%
% IDEA CLAVE: EN OCTAVE UNA SEÑAL DISCRETA SON *DOS* VECTORES
%   Octave no sabe nada de "n". Un vector es solo una lista de números:
%       x = [1 1 1 1 1 1]
%   Para saber en qué instante n está cada valor, guardamos aparte otro
%   vector del mismo largo con los índices de tiempo:
%       nx = [-3 -2 -1 0 1 2]
%   Así, x(1) es el valor en n = nx(1) = -3, x(2) en n = -2, etc.
%   Esto es importante porque conv() devuelve solo los valores, y somos
%   NOSOTROS los que tenemos que calcular a qué n corresponde cada uno.
%
% MINI-GUÍA DE SINTAXIS OCTAVE QUE APARECE EN ESTE ARCHIVO
%   %  ...          -> comentario (Octave ignora el resto de la línea)
%   a = 5;          -> asigna. El ";" final evita que se imprima en consola
%   a:b             -> vector de a hasta b de a 1:  -3:2 = [-3 -2 -1 0 1 2]
%   v(1), v(end)    -> primer y último elemento (¡los índices empiezan en 1!)
%   [a, b]          -> arma un vector juntando a y b
%   v'              -> transpone: convierte un vector fila en columna
% =========================================================================

% clear  -> borra todas las variables de la memoria (arranca limpio)
% clc    -> limpia la ventana de comandos
% close all -> cierra todas las ventanas de gráficos que hubiera abiertas
clear; clc; close all;

% --- Señal de entrada x[n]: pulso rectangular unitario -------------------
% Valores no nulos en un intervalo con n negativos y positivos
%
% "-3:2" genera el vector fila [-3 -2 -1 0 1 2] (6 elementos).
% Elegimos -3..2 porque la consigna pide que incluya n negativos y positivos.
nx = -3:2;                 % intervalo de tiempo discreto de x[n]

% ones(size(nx)) crea un vector de unos del mismo tamaño que nx.
%   size(nx) devuelve [1 6] (1 fila, 6 columnas), así que esto es
%   equivalente a escribir x = [1 1 1 1 1 1].
% Fuera de nx la señal vale 0 implícitamente (no la guardamos).
x  = ones(size(nx));       % x[n] = 1 para -3 <= n <= 2

% --- Respuesta al impulso h[n] = n ----------------------------------------
% nh = [0 1 2 3 4 5]: empieza en 0 y son todos positivos, como pide la
% consigna.
nh = 0:5;                  % comienza en n = 0, todos positivos

% Como h[n] = n, los VALORES de h son iguales a los índices: h = [0 1 2 3 4 5].
% Notar que h[0] = 0, así que el primer valor de h es cero.
h  = nh;                   % h[n] = n

% --- Respuesta del sistema LTI: y[n] = x[n] * h[n] -----------------------
% conv(x, h) calcula la convolución discreta:
%       y[n] = sum_k  x[k] . h[n - k]
% Devuelve un vector de largo  length(x) + length(h) - 1 = 6 + 6 - 1 = 11.
%
% OJO: conv trabaja solo con los valores, "no sabe" dónde empieza cada
% señal en el eje n. Por eso el paso siguiente es construir ny a mano.
y  = conv(x, h);

% REGLA DEL ÍNDICE DE LA CONVOLUCIÓN:
%   si x arranca en n = nx(1) y h arranca en n = nh(1), entonces y arranca
%   en n = nx(1) + nh(1) y termina en n = nx(end) + nh(end).
%   Acá: arranca en -3 + 0 = -3 y termina en 2 + 5 = 7  ->  ny = -3:7
%   (11 elementos, igual que y. ¡Siempre verificar que coincidan!)
% La convolución empieza en nx(1)+nh(1) y termina en nx(end)+nh(end)
ny = (nx(1) + nh(1)) : (nx(end) + nh(end));

% --- Escala común de tiempo discreto: incluye toda la convolución --------
% Usamos el rango de y (el más largo de los tres) para los 3 gráficos,
% con un margen de 1 muestra a cada lado para que no queden pegados al borde.
% lim = [-4 8]
lim = [ny(1) - 1, ny(end) + 1];

% figure abre una ventana de gráficos nueva. 'Name' es el título de la ventana.
figure('Name', 'Ejercicio 1 - Convolucion discreta');

% subplot(filas, columnas, posición) divide la ventana en una grilla y
% elige en qué casillero dibujar. subplot(3,1,1) = 3 filas, 1 columna,
% casillero 1 (arriba de todo).
subplot(3, 1, 1);
% stem(n, v) dibuja una señal discreta: un "palito" vertical con un círculo
% en cada punto (n(i), v(i)). Es el gráfico típico de señales discretas.
% 'filled' rellena los círculos.
stem(nx, x, 'filled');
% xlim fija los límites del eje horizontal. Al usar el mismo lim en los 3
% gráficos, quedan alineados y se pueden comparar visualmente.
% grid on dibuja la cuadrícula de fondo.
xlim(lim); grid on;
title('Entrada x[n] (pulso rectangular)');
ylabel('x[n]');

% Casillero 2 (medio)
subplot(3, 1, 2);
stem(nh, h, 'filled');
xlim(lim); grid on;
title('Respuesta al impulso h[n] = n');
ylabel('h[n]');

% Casillero 3 (abajo)
subplot(3, 1, 3);
stem(ny, y, 'filled');
xlim(lim); grid on;
title('Respuesta y[n] = x[n] * h[n]');
xlabel('n'); ylabel('y[n]');

% Tabla de resultados en consola
% disp muestra algo en la consola.
% ny' e y' son ny e y transpuestos (de fila a columna); [ny' y'] los pone
% lado a lado formando una tabla de 2 columnas: n | y[n].
%
% RESULTADO ESPERADO (para verificar a mano):
%    n : -3 -2 -1  0  1  2  3  4  5  6  7
%  y[n]:  0  1  3  6 10 15 15 14 12  9  5
% Interpretación: mientras la rampa "entra" en el pulso, y[n] crece como
% suma acumulada 0, 0+1, 0+1+2, ...; cuando el pulso cubre toda la rampa
% llega al máximo 0+1+2+3+4+5 = 15 (n = 2 y 3) y luego decrece a medida
% que la rampa "sale" del pulso (se pierden primero los valores chicos).
disp('   n    y[n]');
disp([ny' y']);
