% =========================================================================
% EJERCICIO 4 (continuación): estudios pedidos en el enunciado
%   a) Mejora de H(f) con Hamming frente a la ventana rectangular (Gibbs)
%   b) Banda de transición al modificar el núcleo Nh
%   c) Espectro de la entrada al modificar la ventana no nula Nx
%   d) Calidad de representación al modificar la frecuencia de muestreo Fm
% =========================================================================
clear; clc; close all;

Fc = 500;  Fm = 10000;  Nfft = 8192;      % FFT larga: H(f) y X(f) casi continuos
f  = (0:Nfft-1) * Fm / Nfft;
dB = @(v) 20*log10(abs(v) + 1e-12);

% ---------- a) Rectangular vs Hamming -------------------------------------
Nh = 101;  nh = -(Nh-1)/2:(Nh-1)/2;
h_rect = (2*Fc/Fm) * sinc(2*Fc*nh/Fm);
h_hamm = h_rect .* hamming(Nh)';   h_hamm = h_hamm / sum(h_hamm);
H_rect = fft(h_rect, Nfft);        H_hamm = fft(h_hamm, Nfft);

figure('Name', 'a) Rectangular vs Hamming');
subplot(2, 1, 1);
plot(f, abs(H_rect), f, abs(H_hamm)); xlim([0 Fm/2]); grid on;
legend('Rectangular (Ej. 2 y 3)', 'Hamming normalizada (Ej. 4)');
title(sprintf('|H(f)|, N_h = %d, F_c = %d Hz', Nh, Fc));
subplot(2, 1, 2);
plot(f, dB(H_rect), f, dB(H_hamm)); xlim([0 Fm/2]); ylim([-100 5]); grid on;
title('|H(f)| en dB'); xlabel('f [Hz]'); ylabel('dB');

printf('Rectangular: sobrepico = %.1f %%, lobulo lateral max = %.1f dB\n', ...
       (max(abs(H_rect)) - 1)*100, max(dB(H_rect(f > 1.3*Fc & f < Fm/2))));
printf('Hamming    : sobrepico = %.1f %%, lobulo lateral max = %.1f dB\n', ...
       (max(abs(H_hamm)) - 1)*100, max(dB(H_hamm(f > 1.3*Fc & f < Fm/2))));

% ---------- b) Banda de transición vs Nh -----------------------------------
figure('Name', 'b) Banda de transicion vs Nh');
hold on;
lista_Nh = [21 51 101 201];
for Nh = lista_Nh
  nh = -(Nh-1)/2:(Nh-1)/2;
  h  = (2*Fc/Fm) * sinc(2*Fc*nh/Fm) .* hamming(Nh)';  h = h / sum(h);
  H  = abs(fft(h, Nfft));
  plot(f, H);
  % ancho de transición: de |H| = 0.9 a |H| = 0.1
  f1 = f(find(H < 0.9, 1));  f2 = f(find(H < 0.1, 1));
  printf('Nh = %3d: transicion %4.0f - %4.0f Hz (ancho %4.0f Hz), |H(1000 Hz)| = %.4f\n', ...
         Nh, f1, f2, f2 - f1, H(round(1000*Nfft/Fm) + 1));
end
hold off; xlim([0 2000]); grid on;
legend(arrayfun(@(v) sprintf('N_h = %d', v), lista_Nh, 'UniformOutput', false));
title('|H(f)| con ventana de Hamming para distintos N_h'); xlabel('f [Hz]');

% ---------- c) Espectro de la entrada vs ventana no nula Nx ---------------
figure('Name', 'c) Efecto de enventanado sobre X(f)');
lista_Nx = [101 401 1601];
for i = 1:numel(lista_Nx)
  Nx = lista_Nx(i);  n = 0:Nx-1;
  x  = 10*cos(2*pi*150*n/Fm) + 5*sin(2*pi*1000*n/Fm);
  X  = abs(fft(x, Nfft)) / Nx * 2;       % escalado: el pico tiende a la amplitud
  subplot(numel(lista_Nx), 1, i);
  plot(f, X); hold on;
  stem([150 1000], [10 5], 'r', 'filled'); hold off;   % espectro de x(t) infinita
  xlim([0 1500]); ylim([0 11]); grid on;
  title(sprintf('2|X(f)|/N_x con N_x = %d (ventana de %.1f ms)', Nx, 1000*Nx/Fm));
end
xlabel('f [Hz]');

% ---------- d) Variación de la frecuencia de muestreo Fm ------------------
fd = figure('Name', 'd) Variacion de Fm');
lista_Fm = [10000 3000 2000 1200];
T = 0.04;                                 % misma duración de ventana [s]
for i = 1:numel(lista_Fm)
  Fmi = lista_Fm(i);
  Nx  = round(T*Fmi) + 1;
  Nh  = 2*round(0.005*Fmi) + 1;           % núcleo de ~10 ms en todos los casos
  [t_y, y_conv, ~, fi, X] = ej4(Fmi, Nh, Fc, Nx);  close(gcf);
  figure(fd);                             % volver a la figura del estudio
  t_x = (0:Nx-1) / Fmi;
  x   = 10*cos(2*pi*150*t_x) + 5*sin(2*pi*1000*t_x);

  subplot(numel(lista_Fm), 2, 2*i-1);
  plot(t_x, x, '.-', t_y, y_conv, '.-'); xlim([0 T]); grid on;
  title(sprintf('F_m = %d Hz: x (azul) e y (rojo)', Fmi));
  subplot(numel(lista_Fm), 2, 2*i);
  plot(fi, abs(X)); xlim([0 Fmi/2]); grid on;
  title(sprintf('|X(f)|, F_m = %d Hz (Nyquist = %d Hz)', Fmi, Fmi/2));
end
xlabel('f [Hz]');
