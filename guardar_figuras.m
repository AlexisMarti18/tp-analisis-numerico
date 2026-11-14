% Guarda todas las figuras abiertas como PNG en ./figuras_octave y las cierra
%
% Función auxiliar usada por correr_todo.m. No es parte de la consigna.
% Uso:  guardar_figuras('ej2')  -> genera ej2_1.png, ej2_2.png, ...
function guardar_figuras(prefijo)
  % pwd = carpeta actual. fullfile une partes de una ruta con la barra
  % correcta para el sistema operativo -> '.../TP Laboratorio/figuras_octave'.
  d = fullfile(pwd, 'figuras_octave');
  % Si la carpeta no existe, la crea. "~" es el NO lógico.
  if ~exist(d, 'dir'), mkdir(d); end
  % findall(0, 'type', 'figure') devuelve los identificadores de todas las
  % figuras abiertas; sort los ordena para que se numeren en el orden en
  % que se crearon.
  h = sort(findall(0, 'type', 'figure'));
  for i = 1:numel(h)
    % print guarda la figura h(i) en un archivo:
    %   '-dpng' = formato PNG,  '-r110' = resolución de 110 puntos por pulgada.
    print(h(i), fullfile(d, sprintf('%s_%d.png', prefijo, i)), '-dpng', '-r110');
  end
  close all;   % cierra las figuras para que no se mezclen con el próximo ejercicio
end
