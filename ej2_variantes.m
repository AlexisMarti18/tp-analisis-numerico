% =========================================================================
% EJERCICIO 2 (continuación): se repite el ejercicio cambiando dt y Nh
% =========================================================================
%
% QUÉ PIDE LA CONSIGNA
%   "Resolver el ejercicio nuevamente cambiando el valor de Δt de 0.0001 a
%    0.0005 y a 0.001. También modificando el valor de Nh. Observar las
%    diferencias y sacar conclusiones."
%
% CÓMO SE RESUELVE
%   El caso base (dt = 0.0001, Nh = 101) es el que ya grafica ej2.m, así
%   que acá no se repite. Este script abre solo 4 gráficos:
%     a) dt = 0.0005 y dt = 0.001   (con Nh = 101)     -> 2 gráficos
%     b) Nh = 21 y Nh = 201         (con dt = 0.0001)  -> 2 gráficos
%   Cada llamada a la función ej2(dt, Nh) abre una ventana nueva.
%
% SINTAXIS NUEVA EN ESTE ARCHIVO
%   for v = [a b]        -> repite el bloque con v = a y después con v = b
%     ...
%   end
%   printf(formato, valores) -> imprime texto con formato en la consola
%                               (como sprintf, pero lo muestra directamente).
%       %-7g  -> número, alineado a la izquierda en 7 caracteres
%       %6g   -> número en 6 caracteres
%       %4.2f -> número con 2 decimales
%       \n    -> salto de línea
%   ...      -> al final de una línea: "la instrucción sigue en la próxima"
%   max(abs(y)) -> el mayor valor absoluto de y (la amplitud máxima)
% =========================================================================
clear; clc; close all;

% --- a) Variación del paso de muestreo (Nh = 101, Fc = 500 Hz) -----------
% Para cada dt se muestra en consola:
%   Fm = 1/dt          -> frecuencia de muestreo
%   2*Fc*dt            -> el factor que aparece en h[n] = 2FcΔt sinc(2FcΔt n)
%   max|y|             -> amplitud de la salida
%
% RESULTADOS (corridos en Octave 11.3; el caso base dt = 0.0001 da 10.078):
%   dt = 0.0005  Fm =  2000 Hz  2*Fc*dt = 0.50  max|y| = 10.089
%   dt = 0.001   Fm =  1000 Hz  2*Fc*dt = 1.00  max|y| = 10.000
%
% INTERPRETACIÓN (lo más importante para las conclusiones)
%   * Caso base dt = 0.0001 (Fm = 10000 Hz, en ej2.m): Nyquist = 5000 Hz,
%     ambas componentes bien representadas y el filtro elimina la de
%     1000 Hz. La salida es ~10 cos(2π·150·t).
%
%   * dt = 0.0005 (Fm = 2000 Hz): Nyquist = 1000 Hz, ¡justo la frecuencia de
%     la componente alta! Sus muestras valen 5·sen(2π·1000·n·0.0005) =
%     5·sen(π·n) = 0 para todo n entero. La componente de 1000 Hz
%     DESAPARECE YA AL MUESTREAR, antes de pasar por el filtro.
%
%   * dt = 0.001 (Fm = 1000 Hz): Nyquist = 500 Hz < 1000 Hz -> hay ALIASING.
%     5·sen(2π·1000·n·0.001) = 5·sen(2π·n) = 0: también desaparece.
%     Además 2·Fc·dt = 1, entonces h[n] = sinc(n), que vale 1 en n = 0 y 0
%     en cualquier otro entero: h[n] = δ[n]. El "filtro" es la identidad,
%     y = x (por eso max|y| = 10 exacto). Fc = 500 Hz coincide con Nyquist:
%     el filtro deja pasar todo lo que se puede representar.
%
%   CONCLUSIÓN: en estos dos casos la salida se ve "limpia", pero no
%   gracias al filtro sino al muestreo. Si la componente alta fuera un
%   coseno en vez de un seno, no se anularía: aparecería disfrazada como
%   otra frecuencia (alias), y el filtro ya no podría separarla. Moraleja:
%   hay que muestrear a Fm > 2·(frecuencia máxima de la señal) (Nyquist).
%   Otra diferencia visible: con las mismas 401 muestras, la ventana dura
%   0.04 s, 0.2 s y 0.4 s respectivamente (Nx·dt), así que se ven 6, 30 y
%   60 periodos de la componente de 150 Hz.
for dt = [0.0005 0.001]
  % Llama a ej2 y se queda solo con las dos primeras salidas (ny, y).
  [ny, y] = ej2(dt, 101);
  printf('dt = %-7g  Fm = %6g Hz  2*Fc*dt = %4.2f  max|y| = %6.3f\n', ...
         dt, 1/dt, 2*500*dt, max(abs(y)));
end

% --- b) Variación del largo del núcleo (dt = 0.0001, Fc = 500 Hz) --------
% Se prueba un Nh más chico (21) y uno más grande (201) que el caso base
% (101), para ver el efecto en las dos direcciones.
% Acá no se piden salidas: ej2(...) solo abre los gráficos.
%
% QUÉ OBSERVAR
%   * La sinc con Fc = 500 Hz y dt = 0.0001 cruza por cero cada 10 muestras.
%     Con Nh = 21 solo entra el lóbulo central y un poco más: es una
%     aproximación MUY pobre del filtro ideal. Medido: deja pasar ~6 % de
%     la componente de 1000 Hz (queda una ondulación de ~0.3 de amplitud)
%     y amplifica ~9 % la de 150 Hz (por el sobrepico de Gibbs). En la
%     zona central, la salida se aparta hasta 1.2 de la ideal 10cos(2π150t).
%   * Con Nh = 201 la sinc está mucho más completa: filtro más parecido al
%     ideal (corte más abrupto) y salida más limpia.
%   * Costo de un Nh grande: más cuentas (cada y[n] suma Nh productos) y
%     transitorios más largos en los bordes de y[n] (duran ~Nh muestras),
%     y la salida total es más larga (Nx + Nh - 1).
for Nh = [21 201]
  ej2(0.0001, Nh);
end

% NOTA: el caso Fc = 1000 Hz (el ejemplo del enunciado) no se grafica acá
% porque la consigna no pide variar Fc. Si lo piden en la defensa:
%   >> ej2([], [], 1000)
% La componente de 1000 Hz cae justo en el corte, donde la respuesta en
% frecuencia vale ~0.5, y pasa atenuada a la mitad (amplitud ~2.5).
