% Ejecuta todos los ejercicios y guarda las figuras en ./figuras_octave
% Uso (desde esta carpeta):  octave-cli correr_todo.m   o   >> correr_todo
%
% CÓMO USARLO SI NUNCA USASTE OCTAVE
%   1) Abrir Octave (la versión con interfaz gráfica, "Octave GUI").
%   2) En el panel "Explorador de archivos" (o con el comando cd) ir a la
%      carpeta de este TP. Octave solo encuentra los .m de la carpeta actual.
%   3) En la "Ventana de comandos" escribir el nombre de un archivo SIN el
%      .m y apretar Enter. Por ejemplo:
%         >> ej1            % corre el ejercicio 1 y muestra el gráfico
%         >> ej2            % ejercicio 2 con los valores por defecto
%         >> ej2(0.0005)    % ejercicio 2 con dt = 0.0005
%         >> correr_todo    % corre todo y guarda los PNG (sin mostrarlos)
%   4) Para ver una variable, escribir su nombre sin ";"  (ej.: >> nx).
%      Para pedir ayuda de un comando:  >> help conv
%
% NOTA: este script NO muestra las ventanas (las hace invisibles para
% guardarlas como imagen). Para VER los gráficos, correr cada ejercicio
% por separado (ej1, ej2, ...).

% set(0, propiedad, valor) cambia una configuración por defecto de Octave
% (el 0 significa "la raíz", es decir, valen para todas las figuras):
%   - defaultfigurevisible 'off': las figuras se crean ocultas.
%   - defaultfigureposition [x y ancho alto]: tamaño en píxeles de cada figura.
set(0, 'defaultfigurevisible', 'off');
set(0, 'defaultfigureposition', [50 50 1100 850]);

% Cada línea: muestra un título en consola, corre el ejercicio y guarda
% sus figuras con un prefijo (ej1_1.png, ej2_1.png, ...).
disp('===== EJERCICIO 1 =====');            ej1;            guardar_figuras('ej1');
disp('===== EJERCICIO 2 =====');            ej2;            guardar_figuras('ej2');
disp('===== EJERCICIO 2 - variantes =====');ej2_variantes;  guardar_figuras('ej2_variantes');
disp('===== EJERCICIO 3 =====');            ej3;            guardar_figuras('ej3');
disp('===== EJERCICIO 4 =====');            ej4;            guardar_figuras('ej4');
disp('===== EJERCICIO 4 - estudios =====');  ej4_estudios;   guardar_figuras('ej4_estudios');
disp('Listo: figuras en ./figuras_octave');
