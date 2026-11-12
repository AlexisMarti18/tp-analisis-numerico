% =========================================================================
% ANÁLISIS NUMÉRICO - TRABAJOS DE LABORATORIO
% EJERCICIO 1: Convolución discreta de un pulso rectangular con una rampa
% =========================================================================
clear; clc; close all;

% --- Señal de entrada x[n]: pulso rectangular unitario -------------------
% Valores no nulos en un intervalo con n negativos y positivos
nx = -3:2;                 % intervalo de tiempo discreto de x[n]
x  = ones(size(nx));       % x[n] = 1 para -3 <= n <= 2

% --- Respuesta al impulso h[n] = n ----------------------------------------
nh = 0:5;                  % comienza en n = 0, todos positivos
h  = nh;                   % h[n] = n

% --- Respuesta del sistema LTI: y[n] = x[n] * h[n] -----------------------
y  = conv(x, h);
% La convolución empieza en nx(1)+nh(1) y termina en nx(end)+nh(end)
ny = (nx(1) + nh(1)) : (nx(end) + nh(end));

% --- Escala común de tiempo discreto: incluye toda la convolución --------
lim = [ny(1) - 1, ny(end) + 1];

figure('Name', 'Ejercicio 1 - Convolucion discreta');

subplot(3, 1, 1);
stem(nx, x, 'filled');
xlim(lim); grid on;
title('Entrada x[n] (pulso rectangular)');
ylabel('x[n]');

subplot(3, 1, 2);
stem(nh, h, 'filled');
xlim(lim); grid on;
title('Respuesta al impulso h[n] = n');
ylabel('h[n]');

subplot(3, 1, 3);
stem(ny, y, 'filled');
xlim(lim); grid on;
title('Respuesta y[n] = x[n] * h[n]');
xlabel('n'); ylabel('y[n]');

% Tabla de resultados en consola
disp('   n    y[n]');
disp([ny' y']);
