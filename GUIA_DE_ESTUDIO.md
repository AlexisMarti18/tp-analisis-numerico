# Guía de estudio: TP de Laboratorio de Análisis Numérico

Esta guía explica cada consigna del TP, cómo se resolvió en Octave, qué resultados da y cómo interpretarlos. Está pensada para leerla junto con los archivos `.m`, que tienen comentarios línea por línea.

**Orden de lectura recomendado:** sección 1 (Octave mínimo) → sección 2 (teoría) → un ejercicio a la vez, con su `.m` abierto al lado.

## Índice

1. [Octave en 5 minutos](#1-octave-en-5-minutos)
2. [Teoría necesaria](#2-teoría-necesaria)
3. [Ejercicio 1: convolución discreta](#3-ejercicio-1-convolución-discreta)
4. [Ejercicio 2: filtrado con una sinc](#4-ejercicio-2-filtrado-con-una-sinc)
5. [Ejercicio 3: convolución por FFT](#5-ejercicio-3-convolución-por-fft)
6. [Ejercicio 4: filtro FIR con ventana de Hamming](#6-ejercicio-4-filtro-fir-con-ventana-de-hamming)
7. [Diferencias con el enunciado](#7-diferencias-con-el-enunciado)
8. [Preguntas para preparar la defensa](#8-preguntas-para-preparar-la-defensa)
9. [Mapa de archivos](#9-mapa-de-archivos)

---

## 1. Octave en 5 minutos

### Cómo correr el TP

1. Abrir **Octave (GUI)**.
2. En el explorador de archivos de Octave, ir a la carpeta del TP (Octave solo encuentra los `.m` de la carpeta actual).
3. En la ventana de comandos, escribir el nombre del archivo **sin** `.m` y apretar Enter:

```octave
ej1              % ejercicio 1
ej2              % ejercicio 2 con valores por defecto
ej2(0.0005)      % ejercicio 2 con dt = 0.0005
ej2(0.0001, 21)  % dt = 0.0001 y Nh = 21
ej2_variantes    % todas las variantes que pide la consigna
ej3
ej4
ej4_estudios
correr_todo      % corre todo y guarda los PNG en figuras_octave/
```

### Sintaxis que aparece en el TP

| Escribís | Significa |
|---|---|
| `% texto` | Comentario: Octave lo ignora |
| `a = 5;` | Asignar. El `;` evita que se imprima el resultado |
| `-3:2` | Vector `[-3 -2 -1 0 1 2]` |
| `0:N-1` | Vector `[0 1 2 ... N-1]` |
| `v(1)`, `v(end)` | Primer y último elemento. **Los índices empiezan en 1**, no en 0 |
| `v(1:10)` | Elementos 1 a 10 |
| `[a, b]` | Pegar vectores uno al lado del otro |
| `v'` | Transponer (fila ↔ columna) |
| `ones(1,N)`, `zeros(1,N)` | Vector de N unos / N ceros |
| `a .* b` | Multiplicación **elemento a elemento** (sin el punto es producto matricial) |
| `2 * v` | Multiplica cada elemento de `v` por 2 |
| `sin(v)`, `cos(v)` | Se aplican a cada elemento del vector |
| `sum(v)`, `max(v)`, `abs(v)` | Suma, máximo, módulo/valor absoluto |
| `numel(v)` | Cantidad de elementos |
| `for i = [1 2 3] ... end` | Repetir un bloque |
| `if cond ... end` | Condicional |
| `function [s1,s2] = f(e1,e2) ... end` | Definir una función (el archivo se tiene que llamar `f.m`) |
| `nargin` | Dentro de una función: cuántos parámetros se pasaron |
| `printf('x = %g\n', x)` | Imprimir con formato (`%g` número, `%d` entero, `\n` salto de línea) |
| `@(v) 20*log10(v)` | Función anónima (de una línea) |

### Comandos de señales y gráficos

| Comando | Qué hace |
|---|---|
| `conv(x, h)` | Convolución discreta. Devuelve un vector de largo `length(x)+length(h)-1` |
| `sinc(u)` | `sin(πu)/(πu)`, con `sinc(0) = 1` |
| `fft(x)`, `fft(x, N)` | Transformada rápida de Fourier (con `N`: completa con ceros hasta largo N) |
| `ifft(X)` | FFT inversa |
| `hamming(N)` | Ventana de Hamming de N puntos (como **columna**) |
| `figure` | Abre una ventana de gráficos |
| `subplot(f, c, p)` | Divide la ventana en f×c y dibuja en la posición p (contadas por filas) |
| `stem(n, v)` | Gráfico de señal discreta (palitos) |
| `plot(t, v)` | Gráfico de línea |
| `xlim([a b])` | Límites del eje horizontal |
| `title`, `xlabel`, `ylabel`, `grid on`, `legend` | Decoración del gráfico |

### La idea más importante

**Para Octave una señal es solo una lista de números; no sabe en qué instante n está cada valor.** Por eso cada señal se guarda con dos vectores: valores (`x`) e índices (`nx`). Cuando `conv` devuelve `y`, somos nosotros los que calculamos `ny`.

---

## 2. Teoría necesaria

### 2.1 Señales discretas y sistemas LTI

- Una **señal discreta** x[n] es una secuencia de valores en instantes enteros n.
- Un sistema **LTI** (lineal e invariante en el tiempo) queda completamente caracterizado por su **respuesta al impulso** h[n]: la salida que da cuando la entrada es δ[n] (un 1 en n = 0 y 0 en el resto).
- Para cualquier entrada, la salida es la **convolución**:

  **y[n] = x[n] * h[n] = Σₖ x[k] · h[n − k]**

  Intuitivamente: se da vuelta h, se la desliza sobre x, y en cada posición se multiplica y se suma.

**Reglas de la convolución discreta (se usan en todo el TP):**

- Largo de y = largo de x + largo de h − 1.
- Inicio de y = inicio de x + inicio de h. Fin de y = fin de x + fin de h.

### 2.2 Muestreo, Nyquist y aliasing

- Muestrear x(t) es tomar x[n] = x(n·Δt). **Fm = 1/Δt** es la frecuencia de muestreo.
- **Teorema de Nyquist:** solo se representan correctamente las frecuencias menores a **Fm/2** (frecuencia de Nyquist).
- **Aliasing:** una componente de frecuencia f > Fm/2 aparece "disfrazada" en otra frecuencia, |f − k·Fm|. Ejemplo del TP: 1000 Hz muestreado a 1200 Hz aparece en 200 Hz. Una vez muestreada, **ningún filtro digital puede separarla**.
- Caso particular: una senoidal justo en Fm/2 o en un múltiplo de Fm, si es un **seno**, se muestrea siempre en sus cruces por cero y desaparece.

### 2.3 Filtro pasa bajos ideal

- Deja pasar sin cambios las frecuencias por debajo de **Fc** (frecuencia de corte) y elimina las de arriba. Su respuesta en frecuencia es un rectángulo.
- Su respuesta al impulso es una **sinc**: h(t) = 2Fc·sinc(2Fc·t). Muestreando y multiplicando por Δt:

  **h[n] = 2·Fc·Δt · sinc(2·Fc·Δt·n)**  (la fórmula del enunciado)

- Problema: la sinc es **infinita** y **no causal** (tiene valores para n < 0). En la práctica se la **trunca** a Nh valores centrados en 0.

### 2.4 DFT, FFT y la propiedad de convolución

- La **DFT** pasa una señal de N muestras a N valores complejos X[k], que dicen cuánto hay de cada frecuencia. La **FFT** es un algoritmo rápido para calcularla (más rápido si N = 2^r).
- El índice k corresponde a la frecuencia analógica **f_k = k·Fm/N**. La resolución es Fm/N Hz por cada k.
- Para señales reales, |X[k]| es simétrico respecto de N/2: la mitad superior es el espejo de las frecuencias negativas.
- **Propiedad de convolución:** convolucionar en el tiempo = multiplicar en frecuencia.

  **y = x * h  ⟺  Y[k] = X[k] · H[k]**

- **Detalle importante:** la DFT trata a las señales como periódicas de periodo N, así que el producto da una convolución **circular**. Para que sea igual a la lineal (la de `conv`), hay que completar con ceros hasta **N ≥ Nx + Nh − 1**.

### 2.5 Enventanado, fenómeno de Gibbs y ventana de Hamming

- Cortar una señal a un tramo finito es multiplicarla por una **ventana rectangular**. En frecuencia, eso produce **lóbulos laterales**.
- **Al truncar la sinc (el filtro):** |H(f)| deja de ser un escalón perfecto y aparecen oscilaciones cerca del corte con un sobrepico de ~9 %. Es el **fenómeno de Gibbs**. Agrandar Nh angosta las oscilaciones pero **no reduce el sobrepico**.
- **Al truncar la señal de entrada:** cada raya del espectro (una delta) se convierte en un lóbulo de ancho ~2·Fm/Nx más lóbulos laterales. Se llama **fuga espectral** (*leakage*). Cuanto más larga la ventana, más finos los lóbulos.
- **Ventana de Hamming:** w[m] = 0.54 − 0.46·cos(2πm/(Nh−1)). Es una campana suave que apaga gradualmente los extremos. Elimina el sobrepico y baja mucho los lóbulos laterales, **a cambio** de una banda de transición más ancha.

---

## 3. Ejercicio 1: convolución discreta

**Archivo:** [`ej1.m`](ej1.m) · **Figura:** [`figuras_octave/ej1_1.png`](figuras_octave/ej1_1.png)

### Consigna

Definir x[n] (pulso rectangular unitario con n negativos y positivos) y h[n] = n (desde n = 0). Calcular y[n] con `conv` y graficar las tres señales en tres filas (`subplot` + `stem`), con el mismo eje n (`xlim`) que incluya toda la convolución.

### Resolución

| Señal | Índices | Valores |
|---|---|---|
| x[n] | `nx = -3:2` | `[1 1 1 1 1 1]` |
| h[n] | `nh = 0:5` | `[0 1 2 3 4 5]` |
| y[n] | `ny = -3:7` | `conv(x, h)`, 11 valores |

El eje de y sale de la regla de la sección 2.1: empieza en −3 + 0 = −3 y termina en 2 + 5 = 7.

### Resultado (impreso en consola)

```
 n   : -3 -2 -1  0  1  2  3  4  5  6  7
y[n] :  0  1  3  6 10 15 15 14 12  9  5
```

### Cómo explicarlo

- **Mientras la rampa va "entrando" en el pulso,** y[n] es una suma acumulada: 0, 0+1, 0+1+2...
- **Cuando el pulso cubre toda la rampa,** se llega al máximo 0+1+2+3+4+5 = 15.
- **Después la rampa "sale" del pulso** y se van perdiendo primero los valores chicos (15 − 1 = 14, 14 − 2 = 12...).
- **Verificación de un punto a mano:** y[0] = Σ x[k]·h[0−k], con k de −3 a 0 → h[3] + h[2] + h[1] + h[0] = 6. ✓

---

## 4. Ejercicio 2: filtrado con una sinc

**Archivos:** [`ej2.m`](ej2.m) (función) y [`ej2_variantes.m`](ej2_variantes.m) · **Figuras:** `ej2_1.png` y `ej2_variantes_1.png` a `_7.png`

### Consigna

- Muestrear x(t) = 10 cos(2π·150·t) + 5 sen(2π·1000·t) con Nx = 401 muestras, Δt = 0.0001 s, desde t = 0.
- Armar h[n] = 2·Fc·Δt·sinc(2·Fc·Δt·n) centrada, con Nh impar.
- Calcular y[n] con `conv` y graficar como en el Ej. 1, con el eje limitado a los valores no nulos de y.
- Repetir con Δt = 0.0005 y 0.001, y con otros Nh. Sacar conclusiones.

### Resolución

| Parámetro | Valor | Por qué |
|---|---|---|
| Δt | 0.0001 s (Fm = 10 kHz) | Ejemplo del enunciado |
| Nx | 401 | Ejemplo del enunciado → ventana de 40 ms |
| Fc | **500 Hz** | Tiene que estar *entre* 150 y 1000 Hz (ver sección 7) |
| Nh | 101 | Impar → n = −50 ... 50 |

- `nx = 0:400`, `t = nx*dt`, y x se calcula con `cos` y `sin` sobre todo el vector.
- `nh = -50:50` y `h = 2*Fc*dt*sinc(2*Fc*dt*nh)`. Con estos valores, h[0] = 0.1 y h cruza por cero cada 10 muestras.
- `y = conv(x, h)`. El eje va de 0 + (−50) = **−50** a 400 + 50 = **450**: son 501 valores.
- Se escribió como **función** para poder llamarla con distintos parámetros sin copiar código.

### Qué se observa (caso base)

- **En la zona central,** la salida es una cosenoidal limpia de 150 Hz y amplitud ~10: la componente de 1000 Hz fue eliminada.
- **En los bordes hay transitorios** (subida y bajada gradual, de ~Nh muestras de largo), porque ahí h se superpone solo en parte con x.
- **y empieza en n negativo** porque h está centrada (filtro no causal). La ventaja es que la salida no queda retrasada respecto de la entrada.

### Variantes de Δt (Nh = 101, Fc = 500 Hz)

| Δt | Fm | Nyquist | 2·Fc·Δt | max\|y\| | Qué pasa |
|---|---|---|---|---|---|
| 0.0001 | 10 kHz | 5000 Hz | 0.1 | 10.078 | Caso correcto: el filtro elimina los 1000 Hz |
| 0.0005 | 2 kHz | 1000 Hz | 0.5 | 10.089 | 1000 Hz cae **justo en Nyquist**: 5·sen(π·n) = 0, desaparece **al muestrear** |
| 0.001 | 1 kHz | 500 Hz | 1.0 | 10.000 | **Aliasing**: 5·sen(2π·n) = 0. Además h[n] = sinc(n) = δ[n], o sea **y = x** (el filtro no hace nada) |

**Conclusión clave:** en los dos últimos casos la salida se ve limpia, pero **no gracias al filtro sino al muestreo**. Si la componente alta fuera un coseno, no se anularía: aparecería como un alias y el filtro no podría separarla. Hay que muestrear con **Fm > 2·f_max**.

Otra diferencia: con las mismas 401 muestras, la ventana dura 0.04 s, 0.2 s y 0.4 s, así que se ven 6, 30 y 60 periodos de la componente de 150 Hz.

### Variantes de Nh (Δt = 0.0001)

- **Nh = 21:** entra poco más que el lóbulo central de la sinc, y el filtro es muy pobre. Deja pasar ~6 % de la componente de 1000 Hz y amplifica ~9 % la de 150 Hz (sobrepico de Gibbs). La salida se aparta hasta 1.2 de la ideal.
- **Nh = 51, 101, 201:** cada vez más parecido al filtro ideal, con salida más limpia.
- **Costo de un Nh grande:** más cuentas, transitorios más largos y una salida más larga (Nx + Nh − 1).

### Caso Fc = 1000 Hz (el ejemplo del enunciado)

La componente de 1000 Hz cae justo en el corte, donde |H| = 0.498 (medido). Pasa **atenuada a la mitad**: amplitud ~2.5 en vez de 5. Por eso conviene una Fc bien en el medio.

---

## 5. Ejercicio 3: convolución por FFT

**Archivo:** [`ej3.m`](ej3.m) · **Figura:** [`figuras_octave/ej3_1.png`](figuras_octave/ej3_1.png)

### Consigna

- Resolver el Ej. 2 con la propiedad de convolución de la transformada de Fourier, usando FFT.
- Completar x[n] y h[n] con ceros hasta N = 2^r (ej. 1024), suficiente para y[n].
- Figura de 3×2: tiempo [n] a la izquierda, |FFT| vs [k] a la derecha.
- Comparar con el Ej. 2, notar Gibbs y comparar los espectros de entrada y salida.

### Resolución paso a paso

1. **Mismas x y h del Ej. 2.**
2. **Chequeo de N:** se necesita N ≥ Nx + Nh − 1 = 501, y N = 1024 alcanza. Si no alcanzara, la función corta con un error.
3. **Completar con ceros:** `xN = [x, zeros(1, N-Nx)]` y `hN = [h, zeros(1, N-Nh)]`.
4. **Transformar, multiplicar y antitransformar:** `X = fft(xN)`, `H = fft(hN)`, `Y = X .* H`, `yN = real(ifft(Y))`. El `real` descarta una parte imaginaria de ~10⁻¹⁵ que es solo error de redondeo.
5. **Eje n:** el primer valor de `yN` corresponde a n = 0 + (−50) = −50, igual que con `conv`. Se toman los primeros 501 valores.
6. **Comparación con `conv`:** la máxima diferencia es **1.2×10⁻¹⁴**. Son iguales salvo redondeo de la máquina.

### Cómo leer el eje k

Con Fm = 10000 Hz y N = 1024, cada k equivale a ~9.77 Hz:

| Frecuencia | k |
|---|---|
| 150 Hz | ≈ 15 |
| 500 Hz (Fc) | ≈ 51 |
| 1000 Hz | ≈ 102 |
| 5000 Hz (Nyquist) | 512 |

Lo que hay de k = 512 a 1023 es el espejo de 0–512.

### Qué se ve en cada gráfico

| Posición | Gráfico | Qué observar |
|---|---|---|
| 1 | x[n] con ceros | 401 muestras de señal + 623 ceros |
| 2 | \|X[k]\| | Picos en k ≈ 15 (≈ 1922) y k ≈ 102 (≈ 962), más sus espejos. Tienen "faldas" por la ventana finita (fuga espectral) |
| 3 | h[n] con ceros | La sinc, desde n = −50 |
| 4 | \|H[k]\| | Escalón (≈ 1 hasta k ≈ 51, ≈ 0 después) **con oscilaciones y sobrepico del 10 %**: Gibbs |
| 5 | y[n] = IFFT | Igual al Ej. 2 en los primeros 501 valores; después solo ceros (N alcanzó) |
| 6 | \|Y[k]\| | Sobrevive el pico de 150 Hz y desaparece el de 1000 Hz: el filtrado visto en frecuencia |

La altura teórica de los picos de |X| es amplitud·Nx/2 (2005 y 1002). Da algo menos porque 150 y 1000 Hz no caen exactamente sobre un k entero.

### Conclusiones

- La convolución por FFT da **el mismo resultado** que `conv`. Para señales largas es mucho más rápida: ~N·log N operaciones contra ~Nx·Nh.
- Hace falta el **relleno de ceros** para evitar la convolución circular.
- **El sobrepico de Gibbs en |H| viene de truncar la sinc.** Medido con Nh = 101, 201 y 401: 10.0 %, 9.5 % y 9.2 %. No desaparece agrandando Nh, y por eso el Ej. 4 usa otra ventana.
- **Entrada vs salida:** |Y| = |X|·|H|. Las frecuencias donde |H| ≈ 1 pasan y donde |H| ≈ 0 se eliminan.

---

## 6. Ejercicio 4: filtro FIR con ventana de Hamming

**Archivos:** [`ej4.m`](ej4.m) (función) y [`ej4_estudios.m`](ej4_estudios.m) · **Figuras:** `ej4_1.png` y `ej4_estudios_1.png` a `_4.png`

### Consigna

- Implementar un filtro pasa bajos con la misma señal y sistema.
- Usar **Fm** como dato en lugar de Δt.
- Aplicar una **ventana de Hamming** a h[n] y **normalizar**.
- Calcular la salida por `conv` y por FFT, y comparar.
- Figura de 4×2: tiempo **en segundos** a la izquierda y |espectro| **en Hz** a la derecha. Para el espectro de las dos salidas se usa el mismo Y de la FFT.
- Verificar e investigar: la mejora de Gibbs, la banda de transición vs Nh, el espectro de entrada vs Nx (comparado con la señal infinita) y la representación vs Fm.

### Resolución ([`ej4.m`](ej4.m))

| Cambio respecto de los Ej. 2–3 | Código |
|---|---|
| Fm reemplaza a Δt | `t = n/Fm`, `h = (2*Fc/Fm)*sinc(2*Fc*nh/Fm)`, `f = k*Fm/N` |
| Ventana de Hamming | `h = h .* hamming(Nh)'` (el `'` porque `hamming` devuelve una columna) |
| Normalizado | `h = h / sum(h)` → sum(h) = H(0) = 1 (ganancia 1 en continua) |
| FFT con relleno automático | `fft(x, N)` completa con ceros |
| Ejes físicos | Tiempo en segundos (`n/Fm`) y frecuencia en Hz (`k·Fm/N`), mostrando de 0 a Fm/2 |

¿Por qué normalizar? La ventana cambia la suma de h (acá pasa a 1.0024), y esa suma es la ganancia en continua. Dividiendo por ella, las bajas frecuencias pasan con amplitud exacta.

**Comparación conv vs FFT:** la diferencia máxima es 1.07×10⁻¹⁴, o sea que son idénticas.

### Estudio a) Rectangular vs Hamming (Gibbs)

| Núcleo (Nh = 101) | Sobrepico | Lóbulo lateral máx. (f > 650 Hz) |
|---|---|---|
| Rectangular (Ej. 2–3) | 10.0 % | −27.9 dB |
| Hamming normalizada | **0.0 %** | **−41.0 dB** |

Hamming elimina el sobrepico y deja pasar ~4.5 veces menos amplitud en la banda eliminada. **Gibbs reducido: verificado.**

### Estudio b) Banda de transición vs Nh (con Hamming)

Se mide el ancho de la transición entre |H| = 0.9 y |H| = 0.1:

| Nh | Transición | Ancho | \|H(1000 Hz)\| |
|---|---|---|---|
| 21 | 228 – 989 Hz | 760 Hz | 0.0933 |
| 51 | 316 – 681 Hz | 365 Hz | 0.0017 |
| 101 | 409 – 592 Hz | 183 Hz | 0.0007 |
| 201 | 455 – 546 Hz | 90 Hz | 0.0005 |

- **El ancho es inversamente proporcional a Nh:** duplicar Nh lo reduce a la mitad (Nh × ancho ≈ 18000 ≈ 1.8·Fm).
- **Con Nh = 21 el filtro todavía deja pasar 9 % a 1000 Hz.** Con Nh ≥ 51 ya es despreciable.
- **Un núcleo más largo se acerca más al filtro ideal** (cuyo corte es discontinuo), a costa de más cálculo y transitorios más largos.

### Estudio c) Espectro de la entrada vs largo de la ventana Nx

Se grafica 2|X(f)|/Nx (así el pico tiende a la amplitud) junto a las rayas ideales de la señal infinita (10 en 150 Hz y 5 en 1000 Hz):

| Nx | Duración | Ancho del lóbulo (~2·Fm/Nx) |
|---|---|---|
| 101 | 10 ms | ~200 Hz: picos muy desparramados |
| 401 | 40 ms | ~50 Hz |
| 1601 | 160 ms | ~12 Hz: casi las rayas ideales |

**Efecto del enventanado:** truncar la señal convierte cada raya en un lóbulo con lóbulos laterales. Con una ventana más larga hay **mejor resolución en frecuencia** y el resultado se acerca al espectro de la señal infinita.

### Estudio d) Variación de Fm

La duración de la ventana (40 ms) y del núcleo (~10 ms) se mantienen fijas en segundos; solo cambia la densidad de muestras.

| Fm | Nyquist | Resultado |
|---|---|---|
| 10000 Hz | 5000 Hz | Todo bien representado y suave |
| 3000 Hz | 1500 Hz | 1000 Hz se representa, pero con 3 muestras por periodo: en el tiempo se ve quebrada. El espectro y el filtro siguen funcionando |
| 2000 Hz | 1000 Hz | 1000 Hz cae en Nyquist: el seno da 0 en todas las muestras y **desaparece** |
| 1200 Hz | 600 Hz | **Aliasing:** 1000 Hz aparece en 1200 − 1000 = **200 Hz**. Como 200 < Fc, el filtro **lo deja pasar** y la salida queda contaminada (verificado: hay un pico en 200 Hz de la mitad de altura que el de 150 Hz) |

**Conclusión:** Fm tiene que superar el doble de la máxima frecuencia de la señal. Si no, hay que filtrar en analógico *antes* de muestrear (filtro **antialiasing**). Después del muestreo, el daño no se puede deshacer.

---

## 7. Diferencias con el enunciado

| Enunciado | En el TP | Motivo |
|---|---|---|
| "frecuencias de 30 Hz y 1000 Hz" | 150 Hz | La fórmula dice cos(2π·**150**·t). Error de tipeo de la guía |
| Ejemplo Fc = 1000 Hz | Fc = 500 Hz (1000 Hz se muestra aparte) | La consigna pide Fc *entre* las dos frecuencias. Con 1000 Hz, esa componente pasa a la mitad |
| Comando `haming` | `hamming` | Error de tipeo de la guía |
| "su correspondiente normalizado" | `h / sum(h)` | Interpretación habitual: ganancia unitaria en continua |
| Espectros en Hz | Se muestran de 0 a Fm/2 | La otra mitad es espejo (simetría de la FFT de señales reales) |

---

## 8. Preguntas para preparar la defensa

**¿Por qué hay que calcular `ny` a mano?**
Porque `conv` solo devuelve valores; no sabe en qué n empieza cada señal. Inicio de y = inicio de x + inicio de h.

**¿Por qué Nh tiene que ser impar?**
Para tener una muestra en n = 0 y la misma cantidad a cada lado: h queda simétrica y centrada, y la salida no se retrasa.

**¿Qué significa 2·Fc·Δt en la fórmula de h?**
Es la frecuencia de corte normalizada (Fc dividida por Nyquist). Define cada cuántas muestras cruza por cero la sinc (cada 1/(2FcΔt)) y el valor de h[0].

**¿Por qué con Δt = 0.001 el filtro "no hace nada"?**
Porque 2·Fc·Δt = 1 y sinc(n) es 0 en todo entero salvo n = 0, así que h = δ e y = x. Además Fc = 500 Hz coincide con Nyquist.

**¿Por qué desaparece la componente de 1000 Hz con Δt = 0.0005?**
Porque 1000 Hz es justo Fm/2. Las muestras de un seno en Nyquist caen todas en sus cruces por cero.

**¿Por qué hay que rellenar con ceros antes de la FFT?**
La DFT supone señales periódicas, y el producto X·H da una convolución circular. Con N ≥ Nx + Nh − 1 no hay solapamiento y coincide con la lineal.

**¿Por qué N potencia de 2?**
Porque así el algoritmo FFT es más eficiente.

**¿Qué es el fenómeno de Gibbs y de dónde sale?**
Son las oscilaciones y el sobrepico (~9 %) que aparecen en |H| cerca del corte por truncar la sinc con una ventana rectangular. Agrandar Nh no lo elimina.

**¿Qué mejora la ventana de Hamming y qué empeora?**
Elimina el sobrepico y baja los lóbulos laterales (de −28 a −41 dB acá). A cambio, la transición es más ancha.

**¿Para qué se normaliza h?**
Para que sum(h) = H(0) = 1 y las frecuencias bajas pasen sin cambiar su amplitud.

**¿Cómo se pasa de k a Hz?**
f = k·Fm/N.

**¿Por qué el espectro es simétrico?**
Porque las señales son reales: los k > N/2 representan frecuencias negativas.

**¿Qué pasa si agrando Nh?**
La transición se angosta (proporcional a 1/Nh) y el filtro se acerca al ideal. A cambio hay más cálculo, más transitorio y una salida más larga.

**¿Qué pasa si agrando Nx?**
Los lóbulos del espectro de entrada se angostan (mejor resolución) y se acerca al espectro de la señal infinita.

**¿Qué es el aliasing? Dé un ejemplo del TP.**
Una frecuencia mayor a Fm/2 que aparece como otra más baja. Por ejemplo, 1000 Hz a Fm = 1200 Hz aparece en 200 Hz y pasa el filtro.

**¿Cuál es la diferencia entre `.*` y `*` en Octave?**
`.*` multiplica elemento a elemento; `*` es el producto matricial.

**¿Por qué `real(ifft(...))`?**
Porque el redondeo deja una parte imaginaria de ~10⁻¹⁵ que no tiene significado físico.

---

## 9. Mapa de archivos

| Archivo | Tipo | Contenido |
|---|---|---|
| [`ej1.m`](ej1.m) | script | Ejercicio 1 |
| [`ej2.m`](ej2.m) | función `ej2(dt, Nh, Fc, Nx)` | Ejercicio 2 (un caso) |
| [`ej2_variantes.m`](ej2_variantes.m) | script | Las variantes de Δt, Nh y Fc |
| [`ej3.m`](ej3.m) | función `ej3(dt, Nh, Fc, Nx, N)` | Ejercicio 3 |
| [`ej4.m`](ej4.m) | función `ej4(Fm, Nh, Fc, Nx, N)` | Ejercicio 4 |
| [`ej4_estudios.m`](ej4_estudios.m) | script | Estudios a) a d) del Ej. 4 |
| [`correr_todo.m`](correr_todo.m) | script | Corre todo y guarda los PNG |
| [`guardar_figuras.m`](guardar_figuras.m) | función auxiliar | Guarda las figuras abiertas |
| `figuras_octave/` | imágenes | Figuras generadas |

**Figuras de las variantes:**

- `ej2_variantes_1/2/3` = Δt de 0.0001 / 0.0005 / 0.001.
- `ej2_variantes_4/5/6` = Nh de 21 / 51 / 201.
- `ej2_variantes_7` = Fc = 1000 Hz.
- `ej4_estudios_1/2/3/4` = estudios a / b / c / d.

Todos los números de esta guía salen de correr los scripts en Octave 11.3.
