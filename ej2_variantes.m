% =========================================================================
% EJERCICIO 2 (continuación): se repite el ejercicio cambiando dt y Nh
% =========================================================================
clear; clc; close all;

% --- a) Variación del paso de muestreo (Nh = 101, Fc = 500 Hz) -----------
for dt = [0.0001 0.0005 0.001]
  [ny, y] = ej2(dt, 101);
  printf('dt = %-7g  Fm = %6g Hz  2*Fc*dt = %4.2f  max|y| = %6.3f\n', ...
         dt, 1/dt, 2*500*dt, max(abs(y)));
end

% --- b) Variación del largo del núcleo (dt = 0.0001, Fc = 500 Hz) --------
for Nh = [21 51 201]
  ej2(0.0001, Nh);
end

% --- c) Caso del ejemplo del enunciado: Fc = 1000 Hz ---------------------
% Fc coincide con la frecuencia de la componente de 1000 Hz, por lo que
% esa componente queda en el borde de la banda y pasa atenuada a la mitad.
ej2(0.0001, 101, 1000);
