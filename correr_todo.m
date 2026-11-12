% Ejecuta todos los ejercicios y guarda las figuras en ./figuras_octave
% Uso (desde esta carpeta):  octave-cli correr_todo.m   o   >> correr_todo
set(0, 'defaultfigurevisible', 'off');
set(0, 'defaultfigureposition', [50 50 1100 850]);

disp('===== EJERCICIO 1 =====');            ej1;            guardar_figuras('ej1');
disp('===== EJERCICIO 2 =====');            ej2;            guardar_figuras('ej2');
disp('===== EJERCICIO 2 - variantes =====');ej2_variantes;  guardar_figuras('ej2_variantes');
disp('===== EJERCICIO 3 =====');            ej3;            guardar_figuras('ej3');
disp('===== EJERCICIO 4 =====');            ej4;            guardar_figuras('ej4');
disp('===== EJERCICIO 4 - estudios =====');  ej4_estudios;   guardar_figuras('ej4_estudios');
disp('Listo: figuras en ./figuras_octave');
