% Guarda todas las figuras abiertas como PNG en ./figuras_octave y las cierra
function guardar_figuras(prefijo)
  d = fullfile(pwd, 'figuras_octave');
  if ~exist(d, 'dir'), mkdir(d); end
  h = sort(findall(0, 'type', 'figure'));
  for i = 1:numel(h)
    print(h(i), fullfile(d, sprintf('%s_%d.png', prefijo, i)), '-dpng', '-r110');
  end
  close all;
end
