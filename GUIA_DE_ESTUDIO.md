# Guía de estudio: TP de Laboratorio de Análisis Numérico

Esta guía explica cada consigna del TP, cómo se resolvió en Octave, qué muestra cada gráfico y **por qué** sale así. Está pensada para leerla con Octave abierto: corrés un ejercicio, mirás la ventana "Figure N" y buscás acá la explicación de esa misma figura.

**Orden de lectura recomendado:** sección 1 (Octave mínimo) → sección 2 (teoría) → un ejercicio a la vez, con sus gráficos abiertos al lado.

## Índice

1. [Octave en 5 minutos](#1-octave-en-5-minutos)
2. [Teoría necesaria](#2-teoría-necesaria)
3. [Ejercicio 1: convolución discreta](#3-ejercicio-1-convolución-discreta)
4. [Ejercicio 2: filtrado con una sinc](#4-ejercicio-2-filtrado-con-una-sinc)
5. [Ejercicio 3: convolución por FFT](#5-ejercicio-3-convolución-por-fft)
6. [Ejercicio 4: filtro FIR con ventana de Hamming](#6-ejercicio-4-filtro-fir-con-ventana-de-hamming)
7. [Diferencias con el enunciado](#7-diferencias-con-el-enunciado)
8. [Cómo lo toman: guion de defensa ejercicio por ejercicio](#8-cómo-lo-toman-guion-de-defensa-ejercicio-por-ejercicio)
9. [Preguntas para preparar la defensa](#9-preguntas-para-preparar-la-defensa)
10. [Mapa de archivos y figuras](#10-mapa-de-archivos-y-figuras)

---

## 1. Octave en 5 minutos

### Cómo correr el TP

1. Abrir **Octave (GUI)**.
2. En el explorador de archivos de Octave, ir a la carpeta del TP (Octave solo encuentra los `.m` de la carpeta actual).
3. En la ventana de comandos, escribir el nombre del archivo **sin** `.m` y apretar Enter:

```octave
ej1              % ejercicio 1
close all; ej2   % ejercicio 2 con valores por defecto
ej2_variantes    % las variantes de Δt y Nh que pide la consigna
close all; ej3
close all; ej4
ej4_estudios
correr_todo      % corre todo y guarda los PNG en figuras_octave/
```

### Cómo saber qué figura estás mirando

Cada gráfico se abre en una ventana cuyo título dice **"Figure N: ..."** seguido de una descripción. Por ejemplo:

> Figure 1: Ejercicio 2 - dt = 0.0005 s, Nh = 101, Fc = 500 Hz

En esta guía cada figura se nombra con ese mismo número y título.

- `ej1`, `ej2_variantes` y `ej4_estudios` cierran las ventanas anteriores al empezar, así que **siempre numeran desde Figure 1**.
- `ej2`, `ej3` y `ej4` **no** cierran nada, porque son funciones que se llaman varias veces. Si ya había ventanas abiertas, la nueva va a tener otro número. Por eso conviene escribir `close all;` antes, como en los comandos de arriba.

**Tip:** con 400 muestras, los gráficos de `stem` se ven como una "banda" de palitos. Para ver cada muestra, usá la **lupa con +** de la barra de herramientas de la ventana y hacé zoom en una zona.

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

**Otra forma de verlo, muy útil para el Ej. 2:** si h es una "campanita" de valores positivos, cada y[n] es un **promedio ponderado** de las muestras de x cercanas a n. Promediar suaviza: las oscilaciones rápidas se cancelan (los valores positivos y negativos se compensan dentro del promedio) y las lentas sobreviven. Eso es un filtro pasa bajos.

### 2.2 Muestreo, Nyquist y aliasing

- Muestrear x(t) es tomar x[n] = x(n·Δt). **Fm = 1/Δt** es la frecuencia de muestreo.
- **Teorema de Nyquist:** solo se representan correctamente las frecuencias menores a **Fm/2** (frecuencia de Nyquist).
- **Aliasing:** una componente de frecuencia f > Fm/2 aparece "disfrazada" en otra frecuencia, |f − k·Fm|. Ejemplo del TP: 1000 Hz muestreado a Fm = 1200 Hz aparece en 200 Hz. Una vez muestreada, **ningún filtro digital puede separarla**.
- Caso particular: una senoidal justo en Fm/2 o en un múltiplo de Fm, si es un **seno**, se muestrea siempre en sus cruces por cero y desaparece.

### 2.3 Filtro pasa bajos ideal

- Deja pasar sin cambios las frecuencias por debajo de **Fc** (frecuencia de corte) y elimina las de arriba. Su respuesta en frecuencia es un rectángulo.
- Su respuesta al impulso es una **sinc**: h(t) = 2Fc·sinc(2Fc·t). Muestreando y multiplicando por Δt:

  **h[n] = 2·Fc·Δt · sinc(2·Fc·Δt·n)**  (la fórmula del enunciado)

- Propiedades de esta h[n] que se ven en los gráficos:
  - **Valor máximo:** h[0] = 2·Fc·Δt.
  - **Cruces por cero:** cada 1/(2·Fc·Δt) muestras. En segundos, cada 1/(2·Fc), o sea cada 0.001 s con Fc = 500 Hz, sin importar Δt.
  - **Lóbulo central:** del primer cero negativo al primer cero positivo. Ahí están los pesos grandes, todos positivos.
- Problema: la sinc es **infinita** y **no causal** (tiene valores para n < 0). En la práctica se la **trunca** a Nh valores centrados en 0.

### 2.4 DFT, FFT y la propiedad de convolución

- La **DFT** pasa una señal de N muestras a N valores complejos X[k], que dicen cuánto hay de cada frecuencia. La **FFT** es un algoritmo rápido para calcularla (más rápido si N = 2^r).
- El índice k corresponde a la frecuencia analógica **f_k = k·Fm/N**. La resolución es Fm/N Hz por cada k.
- Para señales reales, |X[k]| es simétrico respecto de N/2: la mitad superior es el espejo de las frecuencias negativas.
- Una senoidal de amplitud A que dura Nx muestras da en |X[k]| un pico de altura ≈ **A·Nx/2**.
- **Propiedad de convolución:** convolucionar en el tiempo = multiplicar en frecuencia.

  **y = x * h  ⟺  Y[k] = X[k] · H[k]**

- **Detalle importante:** la DFT trata a las señales como periódicas de periodo N, así que el producto da una convolución **circular**. Para que sea igual a la lineal (la de `conv`), hay que completar con ceros hasta **N ≥ Nx + Nh − 1**.

### 2.5 Enventanado, fenómeno de Gibbs y ventana de Hamming

- Cortar una señal a un tramo finito es multiplicarla por una **ventana rectangular**. En frecuencia, eso produce **lóbulos laterales**.
- **Al truncar la sinc (el filtro):** |H(f)| deja de ser un escalón perfecto. Aparecen ondulaciones en la banda de paso y un sobrepico de ~9 % justo antes del corte. Es el **fenómeno de Gibbs**. Agrandar Nh angosta las ondulaciones pero **no reduce el sobrepico**.
- **Al truncar la señal de entrada:** cada raya del espectro (una delta) se convierte en un lóbulo de ancho ~2·Fm/Nx más lóbulos laterales. Se llama **fuga espectral** (*leakage*). Cuanto más larga la ventana, más finos los lóbulos.
- **Ventana de Hamming:** w[m] = 0.54 − 0.46·cos(2πm/(Nh−1)). Es una campana suave que apaga gradualmente los extremos. Elimina el sobrepico y baja mucho los lóbulos laterales, **a cambio** de una banda de transición más ancha.

### 2.6 Decibeles (para la Figure 1 de `ej4_estudios`)

dB = 20·log₁₀(|H|). Los valores de referencia:

| \|H\| | 1 | 0.1 | 0.01 | 0.001 | 0.0001 |
|---|---|---|---|---|---|
| dB | 0 | −20 | −40 | −60 | −80 |

Cada −20 dB es **10 veces menos amplitud**. Sirve para ver valores chiquitos que en escala lineal parecen 0.

---

## 3. Ejercicio 1: convolución discreta

**Archivo:** [`ej1.m`](ej1.m) · **Se corre con:** `ej1` · **Abre:** Figure 1 · **PNG:** [`figuras_octave/ej1_1.png`](figuras_octave/ej1_1.png)

### Consigna

Definir x[n] (pulso rectangular unitario con n negativos y positivos) y h[n] = n (desde n = 0). Calcular y[n] con `conv` y graficar las tres señales en tres filas (`subplot` + `stem`), con el mismo eje n (`xlim`) que incluya toda la convolución.

### Resolución

| Señal | Índices | Valores |
|---|---|---|
| x[n] | `nx = -3:2` | `[1 1 1 1 1 1]` |
| h[n] | `nh = 0:5` | `[0 1 2 3 4 5]` |
| y[n] | `ny = -3:7` | `conv(x, h)`, 11 valores |

El eje de y sale de la regla de la sección 2.1: empieza en −3 + 0 = −3 y termina en 2 + 5 = 7.

### Cómo leer la Figure 1: "Ejercicio 1 - Convolucion discreta"

Los tres gráficos tienen el mismo eje horizontal, de **−4 a 8**: el rango de y (−3 a 7) más un lugar de margen a cada lado.

**Arriba, "Entrada x[n] (pulso rectangular)":** 6 palitos de altura 1, en n = −3, −2, −1, 0, 1 y 2. En cualquier otro n la señal vale 0, y por eso no se dibuja nada.

**Medio, "Respuesta al impulso h[n] = n":** un punto sobre el eje en n = 0 (porque h[0] = 0) y palitos de altura 1, 2, 3, 4 y 5 en n = 1 a 5. Es una rampa que **termina** en n = 5.

**Abajo, "Respuesta y[n] = x[n] * h[n]":** 11 valores, de n = −3 a n = 7:

```
 n   : -3 -2 -1  0  1  2  3  4  5  6  7
y[n] :  0  1  3  6 10 15 15 14 12  9  5
```

(Estos valores también se imprimen en la consola.)

### Por qué y[n] tiene esa forma

Como x es todo unos, cada y[n] es simplemente **la suma de los valores de la rampa que caen dentro de una "ventana" de 6 lugares**, y esa ventana avanza un lugar por cada n:

```
y[-2] = 1    [ .  .  .  .  0  1 ]          0+1
y[-1] = 3    [ .  .  .  0  1  2 ]          0+1+2
y[ 0] = 6    [ .  .  0  1  2  3 ]
y[ 1] = 10   [ .  0  1  2  3  4 ]
y[ 2] = 15   [ 0  1  2  3  4  5 ]   <- la rampa entera adentro
y[ 3] = 15   [ 1  2  3  4  5  0 ]   sale el 0 (no resta nada)
y[ 4] = 14   [ 2  3  4  5  0  0 ]   sale el 1
y[ 5] = 12   [ 3  4  5  0  0  0 ]   sale el 2
y[ 6] = 9    [ 4  5  0  0  0  0 ]   sale el 3
y[ 7] = 5    [ 5  0  0  0  0  0 ]   sale el 4
```

1. **Sube (n = −2 a 2):** la rampa va entrando y cada paso suma un valor nuevo.
2. **Máximo 15 (n = 2 y 3):** está toda la rampa adentro: 0+1+2+3+4+5. En n = 3 se pierde el 0, que no resta nada, y por eso se repite el 15.
3. **Baja (n = 4 a 7):** la rampa se terminó, así que no entra nada nuevo y solo se pierden valores. Primero salen los chicos, por eso al principio baja despacio (15 → 14 → 12) y después más rápido (9 → 5).

**Por qué y[−3] = 0 aunque x[−3] = 1:** porque h[0] = 0. Como h vale 0 en n ≤ 0, y[n] solo depende de valores de x **anteriores** a n. En n = −3 no hay nada anterior.

**Por qué hay salida en n = 3 a 7, donde x ya vale 0:** porque y[n] no es x[n]·h[n] en el mismo instante. La rampa corrida todavía "alcanza" los unos de x en n = 1 y 2. Es como el eco de una sala: el sonido sigue después de que terminó el aplauso.

**Verificación a mano de un punto:** y[−1] = x[−3]·h[2] + x[−2]·h[1] = 2 + 1 = 3. ✓

---

## 4. Ejercicio 2: filtrado con una sinc

**Archivos:** [`ej2.m`](ej2.m) (función) y [`ej2_variantes.m`](ej2_variantes.m)

| Se corre con | Abre | Qué muestra | PNG |
|---|---|---|---|
| `close all; ej2` | Figure 1 | Caso base: Δt = 0.0001 s, Nh = 101 | `ej2_1.png` |
| `ej2_variantes` | Figure 1 | Δt = 0.0005 s | `ej2_variantes_1.png` |
| | Figure 2 | Δt = 0.001 s | `ej2_variantes_2.png` |
| | Figure 3 | Nh = 21 | `ej2_variantes_3.png` |
| | Figure 4 | Nh = 201 | `ej2_variantes_4.png` |

### Consigna

- Muestrear x(t) = 10 cos(2π·150·t) + 5 sen(2π·1000·t) con Nx = 401 muestras, Δt = 0.0001 s, desde t = 0.
- Armar h[n] = 2·Fc·Δt·sinc(2·Fc·Δt·n) centrada, con Nh impar.
- Calcular y[n] con `conv` y graficar como en el Ej. 1, con el eje limitado a los valores no nulos de y.
- Repetir con Δt = 0.0005 y 0.001, y con otros Nh. Sacar conclusiones.

### Resolución

| Parámetro | Valor | Por qué |
|---|---|---|
| Δt | 0.0001 s (Fm = 10000 Hz) | Ejemplo del enunciado |
| Nx | 401 | Ejemplo del enunciado → dura 401 × 0.0001 = 0.04 s |
| Fc | **500 Hz** | Tiene que estar *entre* 150 y 1000 Hz (ver sección 7) |
| Nh | 101 | Impar → n = −50 ... 50 |

- `nx = 0:400`, `t = nx*dt`, y x se calcula con `cos` y `sin` sobre todo el vector.
- `nh = -50:50` y `h = 2*Fc*dt*sinc(2*Fc*dt*nh)`.
- `y = conv(x, h)`. El eje va de 0 + (−50) = **−50** a 400 + 50 = **450**: son 501 valores.
- Se escribió como **función** para poder llamarla con distintos parámetros sin copiar código.

### Las dos componentes de la entrada

| | Amplitud | Frecuencia | Periodo | Muestras por periodo con Δt = 0.0001 s | Ciclos en 401 muestras |
|---|---|---|---|---|---|
| 10 cos(2π·150·t) | 10 | 150 Hz | 0.00667 s | ≈ 67 | 6 |
| 5 sen(2π·1000·t) | 5 | 1000 Hz | 0.001 s | 10 | 40 |

---

### Figure 1 de `ej2`: "Ejercicio 2 - dt = 0.0001 s, Nh = 101, Fc = 500 Hz" (caso base)

El eje horizontal va de **−50 a 450** en los tres gráficos: son los límites de y.

**Arriba, "Entrada x[n] (N_x = 401, Δt = 0.0001 s, F_m = 10000 Hz)":**
- Se ven **6 ondas grandes**: el coseno de 150 Hz, de amplitud 10.
- Montada encima hay una **ondulación rápida y chica** de 10 muestras por ciclo: el seno de 1000 Hz, de amplitud 5.
- Los picos llegan a **≈ 14.75**, cerca de 15, porque a veces los máximos de las dos componentes coinciden (10 + 5).
- **x[0] = 10.** En t = 0 el coseno vale su máximo y el seno vale 0. Así se distingue el cos del sen.
- x solo existe de n = 0 a 400. Antes y después no hay nada dibujado.

**Medio, "Respuesta al impulso h[n] = 2F_cΔt sinc(2F_cΔt n) (N_h = 101, F_c = 500 Hz)":**
- **Máximo 0.1 en n = 0** (2·500·0.0001 = 0.1).
- **Cruza por cero cada 10 muestras** (1/(2·500·0.0001) = 10). El lóbulo central va de n = −10 a 10.
- **A los costados hay lóbulos alternados** (−0.02, +0.013, ...) que se van achicando.
- **En n = ±50 se corta de golpe:** ahí terminan las 101 muestras. Ese corte es el que después produce el fenómeno de Gibbs (Ej. 3).

**Abajo, "Respuesta y[n] = x[n] * h[n]":**
- **De n ≈ 50 a 350:** una cosenoidal **suave** de 150 Hz y amplitud 10. La ondulación rápida de 1000 Hz **desapareció**.
- **Primer y último pico más bajos** (≈ 9 y ≈ 8), y pequeñas ondulaciones de n = −50 a 0 y de 400 a 450. Son los **transitorios**.

**Por qué desaparecen los 1000 Hz y quedan los 150 Hz:**
- El lóbulo central de h abarca 20 muestras (n = −10 a 10), y cada y[n] es casi un **promedio ponderado de 20 muestras vecinas de x**.
- **La componente de 1000 Hz completa 2 ciclos en esas 20 muestras.** Sus partes positivas y negativas se compensan y el promedio da ≈ 0.
- **La de 150 Hz tarda 67 muestras en un ciclo.** En 20 muestras casi no cambia, así que el promedio la deja pasar entera.
- En términos de frecuencia: 150 Hz < Fc = 500 Hz pasa, y 1000 Hz > Fc se elimina.

**Por qué hay transitorios:** en los bordes, h se superpone con x solo en parte (fuera de 0–400, x vale 0). El promedio incluye ceros y la salida sube o baja gradualmente.

**Por qué y empieza en n = −50:** h está centrada (tiene valores en n negativos), y por la regla del inicio, 0 + (−50) = −50. La ventaja es que la salida no queda retrasada: los picos de y están en los mismos n que los de la onda de 150 Hz en x.

---

### Figure 1 de `ej2_variantes`: "Ejercicio 2 - dt = 0.0005 s, Nh = 101, Fc = 500 Hz"

Cálculos previos: **Fm = 1/0.0005 = 2000 Hz**, **Nyquist = 1000 Hz**, **2·Fc·Δt = 0.5**.

**Arriba, "Entrada x[n] (N_x = 401, Δt = 0.0005 s, F_m = 2000 Hz)":**
- Parece una banda densa de palitos que van de **−10 a 10**: son 30 ciclos del coseno de 150 Hz (las 401 muestras ahora duran 0.2 s) con ≈ 13 muestras por ciclo. Con zoom se ve la cosenoidal.
- **Los picos llegan a 10, no a 15 como en el caso base.** Eso es la prueba de que **la componente de 1000 Hz no está**.

**Por qué no está la componente de 1000 Hz:** con Δt = 0.0005 s, sus muestras valen 5·sen(2π·1000·n·0.0005) = 5·sen(π·n). sen(π·n) es **0 para todo n entero**: el muestreo cae justo en los cruces por cero del seno. 1000 Hz es exactamente la frecuencia de Nyquist (Fm/2 = 1000 Hz). La componente **desaparece al muestrear, antes de llegar al filtro**.

**Medio, h[n] con N_h = 101:**
- **Máximo 0.5 en n = 0** (2·500·0.0005 = 0.5).
- **h[±1] = 0.32, h[±2] = 0, h[±3] = −0.11.** Cruza por cero **cada 2 muestras** (1/0.5 = 2).
- Por eso se ve como un pico angosto con colitas que se apagan rápido, hasta n = ±50.
- **Por qué ahora es tan angosta:** en **segundos** es la misma sinc del caso base (cruza por cero cada 0.001 s). Pero ahora cada muestra dura 0.0005 s en vez de 0.0001 s, así que 0.001 s son solo 2 muestras.

**Abajo, y[n]:**
- De n = 0 a 400, **igual a x** (diferencia máxima 0.001): pasa el coseno de 150 Hz de amplitud 10.
- Fuera de 0–400 hay un transitorio chico, de hasta ≈ 2.6 justo al lado de los bordes, que se apaga enseguida.

**Conclusión de esta figura:** la salida queda "limpia", pero **no gracias al filtro**: la componente de 1000 Hz ya se había perdido al muestrear.

---

### Figure 2 de `ej2_variantes`: "Ejercicio 2 - dt = 0.001 s, Nh = 101, Fc = 500 Hz"

Cálculos previos: **Fm = 1000 Hz**, **Nyquist = 500 Hz**, **2·Fc·Δt = 1**.

**Arriba, "Entrada x[n] (N_x = 401, Δt = 0.001 s, F_m = 1000 Hz)":**
- Otra vez una banda de **−10 a 10**: ahora son 60 ciclos del coseno de 150 Hz (las 401 muestras duran 0.4 s) con solo ≈ 6.7 muestras por ciclo.
- **Picos de 10:** la componente de 1000 Hz tampoco está.

**Por qué no está:** 1000 Hz > Nyquist = 500 Hz, así que hay **aliasing**. Sus muestras valen 5·sen(2π·1000·n·0.001) = 5·sen(2π·n) = 0 para todo n: el seno da un ciclo completo exacto entre muestra y muestra y siempre se lo "agarra" en 0. Si fuera un coseno, valdría 5 en todas las muestras y aparecería como una constante.

**Medio, h[n] con N_h = 101:**
- **Un único palito de altura 1 en n = 0.** El resto es exactamente 0: la línea gruesa sobre el eje son 100 ceros.

**Por qué h es un solo palito:** h[n] = 1·sinc(n), y sinc(n) = sen(π·n)/(π·n) vale 1 en n = 0 y **0 en cualquier otro entero**. Entonces h[n] = δ[n], el impulso. Convolucionar con δ deja la señal igual, así que este "filtro" no hace nada.

**Abajo, y[n]:**
- De n = 0 a 400, **idéntica a x** (diferencia de 10⁻¹⁵, solo redondeo).
- De n = −50 a −1 y de 401 a 450, **exactamente 0** (las líneas planas sobre el eje). No hay transitorio porque h no tiene colas.

**Conclusión de esta figura:** Fc = 500 Hz coincide con Nyquist, así que el filtro deja pasar todo lo que se puede representar. La componente de 1000 Hz no está por culpa del muestreo.

---

### Figure 3 de `ej2_variantes`: "Ejercicio 2 - dt = 0.0001 s, Nh = 21, Fc = 500 Hz"

Mismo muestreo que el caso base, pero con un núcleo de solo 21 valores.

**Arriba, x[n]:** la misma entrada que en el caso base, con las 6 ondas de 150 Hz, la ondulación de 1000 Hz y picos de ≈ 15.

**Medio, "h[n] ... (N_h = 21, F_c = 500 Hz)":**
- Solo n = −10 a 10: **es únicamente el lóbulo central de la sinc**, todo positivo. Máximo 0.1 en n = 0, que baja a 0 justo en n = ±10.
- **Faltan todos los lóbulos laterales** que se veían en el caso base.

**Abajo, y[n]:**
- La cosenoidal de 150 Hz, pero con amplitud **≈ 10.9** en vez de 10: los picos pasan la línea de 10.
- **Tiene "bultitos"** en las crestas y los valles: queda un resto de la componente de 1000 Hz, de amplitud ≈ 0.29.
- **Los transitorios son cortos:** solo ≈ 10 muestras en cada borde (y existe de n = −10 a 410).

**Por qué filtra peor:**
- Sin los lóbulos laterales, h ya no se parece a la sinc ideal, y su respuesta en frecuencia no es un escalón limpio.
- **A 150 Hz la ganancia es 1.09,** y por eso la salida queda en 10.9.
- **A 1000 Hz la ganancia no llega a 0:** deja pasar ≈ 6 %, o sea 0.29 de los 5 originales.

Intuitivamente: promediar con una ventanita de solo 20 muestras no cancela del todo la componente de 1000 Hz.

**Por qué los transitorios son cortos:** el núcleo dura 21 muestras, así que el efecto de borde dura ≈ 10 muestras de cada lado.

---

### Figure 4 de `ej2_variantes`: "Ejercicio 2 - dt = 0.0001 s, Nh = 201, Fc = 500 Hz"

**Arriba, x[n]:** la misma entrada. El eje ahora va de **−100 a 500**, porque y es más larga.

**Medio, "h[n] ... (N_h = 201, F_c = 500 Hz)":**
- La sinc de n = −100 a 100. El mismo lóbulo central (máximo 0.1, ceros cada 10 muestras), pero ahora con **10 lóbulos laterales de cada lado** en vez de 5.
- Es una aproximación más completa de la sinc ideal.

**Abajo, y[n]:**
- La cosenoidal de 150 Hz **muy limpia**: el resto de 1000 Hz es ≈ 0.035, la mitad que en el caso base.
- **Amplitud ≈ 10.2** (los picos pasan apenas la línea de 10).
- **Transitorios largos:** pequeñas ondulaciones de n = −100 a 0 y de 400 a 500.

**Por qué la amplitud da 10.2 y no exactamente 10:** aunque el núcleo sea más largo, sigue estando **cortado de golpe** (ventana rectangular). Eso deja ondulaciones de ±5 % en la ganancia de la banda que pasa: la curva azul de la Figure 1 de `ej4_estudios` sube y baja entre 0.95 y 1.05. Según dónde caigan los 150 Hz sobre esa ondulación, la ganancia queda un poco arriba o abajo de 1: con Nh = 101 cae en 1.002 y con Nh = 201 en 1.022.

**Conclusión: agrandar Nh mejora el rechazo de los 1000 Hz y hace el corte más abrupto, pero no elimina las ondulaciones de Gibbs.** Para eso hace falta otra ventana, que es lo que hace el Ej. 4.

**Costo de un Nh grande:** más cuentas (cada y[n] suma 201 productos), transitorios más largos y una salida más larga (Nx + Nh − 1 = 601 valores).

---

### Resumen de las variantes

| Figura | Δt | Fm | Nyquist | Nh | Amplitud de y (150 Hz) | Resto de 1000 Hz | Por qué |
|---|---|---|---|---|---|---|---|
| Figure 1 de `ej2` | 0.0001 s | 10000 Hz | 5000 Hz | 101 | 10.0 | 0.07 | El filtro funciona |
| Figure 1 de `ej2_variantes` | 0.0005 s | 2000 Hz | 1000 Hz | 101 | 10.0 | 0 | El seno de 1000 Hz cae en Nyquist y desaparece al muestrear |
| Figure 2 de `ej2_variantes` | 0.001 s | 1000 Hz | 500 Hz | 101 | 10.0 | 0 | Aliasing (desaparece) y h = δ (y = x) |
| Figure 3 de `ej2_variantes` | 0.0001 s | 10000 Hz | 5000 Hz | 21 | 10.9 | 0.29 | Núcleo muy corto: filtro pobre |
| Figure 4 de `ej2_variantes` | 0.0001 s | 10000 Hz | 5000 Hz | 201 | 10.2 | 0.035 | Núcleo largo: rechaza mejor, pero Gibbs sigue |

**Conclusiones para presentar:**

1. **Δt (Fm) tiene que cumplir Nyquist:** Fm > 2 × 1000 = 2000 Hz. Si no, la componente alta desaparece o aparece como otra frecuencia (alias), y el filtro no la puede tratar.
2. **Con las mismas 401 muestras, un Δt más grande abarca más tiempo** (0.04 s, 0.2 s y 0.4 s). Se ven más ciclos y con menos muestras por ciclo.
3. **Un Nh más grande** da un filtro más parecido al ideal: rechaza mejor y corta más abrupto. A cambio, más cálculo y transitorios más largos. El sobrepico de Gibbs no desaparece.

### Caso Fc = 1000 Hz (el ejemplo del enunciado)

La consigna no pide variar Fc, así que este caso no se grafica en `ej2_variantes.m`. Si lo piden en la defensa, se corre con `close all; ej2([], [], 1000)`.

La componente de 1000 Hz cae justo en el corte, donde la ganancia vale 0.498 (medido), y pasa **atenuada a la mitad**: amplitud ≈ 2.5 en vez de 5. En el gráfico de y se ve la cosenoidal de 150 Hz con una ondulación rápida chica encima. Por eso conviene una Fc bien en el medio entre 150 y 1000 Hz.

---

## 5. Ejercicio 3: convolución por FFT

**Archivo:** [`ej3.m`](ej3.m) · **Se corre con:** `close all; ej3` · **Abre:** Figure 1 · **PNG:** [`figuras_octave/ej3_1.png`](figuras_octave/ej3_1.png)

### Consigna

- Resolver el Ej. 2 con la propiedad de convolución de la transformada de Fourier, usando FFT.
- Completar x[n] y h[n] con ceros hasta N = 2^r (ej. 1024), suficiente para y[n].
- Figura de 3×2: tiempo [n] a la izquierda, |FFT| vs [k] a la derecha.
- Comparar con el Ej. 2, notar Gibbs y comparar los espectros de entrada y salida.

### Resolución paso a paso

1. **Mismas x y h del caso base del Ej. 2** (Δt = 0.0001 s, Nh = 101, Fc = 500 Hz).
2. **Chequeo de N:** se necesita N ≥ Nx + Nh − 1 = 501, y N = 1024 alcanza. Si no alcanzara, la función corta con un error.
3. **Completar con ceros:** `xN = [x, zeros(1, N-Nx)]` y `hN = [h, zeros(1, N-Nh)]`.
4. **Transformar, multiplicar y antitransformar:** `X = fft(xN)`, `H = fft(hN)`, `Y = X .* H`, `yN = real(ifft(Y))`. El `real` descarta una parte imaginaria de ~10⁻¹⁵ que es solo error de redondeo.
5. **Eje n:** el primer valor de `yN` corresponde a n = 0 + (−50) = −50, igual que con `conv`. Se toman los primeros 501 valores.
6. **Comparación con `conv`:** la consola imprime `Max |y_fft - y_conv| = 1.243e-14`. Son iguales salvo redondeo de la máquina.

### Cómo leer el eje k

Con Fm = 10000 Hz y N = 1024, cada k equivale a 10000/1024 ≈ 9.77 Hz:

| Frecuencia | k | Espejo (1024 − k) |
|---|---|---|
| 150 Hz | ≈ 15 | ≈ 1009 |
| 500 Hz (Fc) | ≈ 51 | ≈ 973 |
| 1000 Hz | ≈ 102 | ≈ 922 |
| 5000 Hz (Nyquist) | 512 | (centro) |

Todo lo que hay de k = 512 a 1023 es el **espejo** de 0–512, porque representa las frecuencias negativas. Por eso los espectros de esta figura son simétricos.

### Cómo leer la Figure 1: "Ejercicio 3 - Convolucion por FFT (N = 1024)"

Los gráficos de la izquierda tienen el eje n de **−50 a 973** (los 1024 lugares). Los de la derecha tienen el eje k de **0 a 1023**.

**Fila 1, izquierda, "x[n] completada con ceros":** la misma entrada del Ej. 2 de n = 0 a 400 (6 ondas con la ondulación de 1000 Hz, picos de ≈ 15) y después **una línea plana en 0 hasta el final**: son los 623 ceros agregados.

**Fila 1, derecha, "|X[k]|":**
- **Pico alto de ≈ 1922 en k = 15:** el coseno de 150 Hz.
- **Pico de ≈ 962 en k = 102:** el seno de 1000 Hz.
- **Las mismas alturas en k = 1009 y k = 922:** los espejos.
- **Cerca de cada pico, una "falda"** de valores chicos que no llega a 0.

**Por qué esas alturas:**
- **Una senoidal de amplitud A en Nx muestras da un pico de ≈ A·Nx/2.** El de 150 Hz debería medir 10·401/2 = 2005 y el de 1000 Hz, 5·401/2 = 1002.
- **El de 150 Hz es el doble del de 1000 Hz** porque su amplitud es el doble (10 contra 5).
- **Dan algo menos (1922 y 962)** porque 150 Hz y 1000 Hz no caen exactamente sobre un k entero: 150/9.77 = 15.36 y 1000/9.77 = 102.4. Parte de la energía se reparte entre los k vecinos.

**Por qué hay faldas:** la señal está cortada a 401 muestras (ventana rectangular). En vez de rayas perfectas, cada pico queda con lóbulos laterales: es la **fuga espectral** (sección 2.5).

**Fila 2, izquierda, "h[n] completada con ceros":** la sinc del Ej. 2 (máximo 0.1) aparece comprimida a la izquierda, de n = −50 a 50, y después vienen 923 ceros.

**Fila 2, derecha, "|H[k]| (fenómeno de Gibbs en el corte)":**
- **De k = 0 a ≈ 51:** vale ≈ 1, pero **ondulando** (sube y baja alrededor de 1).
- **Justo antes del corte, en k ≈ 41 (≈ 400 Hz), hay un sobrepico de ≈ 1.1.**
- **Después de k ≈ 51, cae a ≈ 0** y se queda ahí.
- **Lo mismo en espejo,** de k ≈ 973 a 1023.

**Por qué |H| tiene forma de escalón:** H es la respuesta en frecuencia del filtro. Vale 1 donde deja pasar (por debajo de Fc = 500 Hz, k ≈ 51) y 0 donde elimina.

**Por qué ondula y tiene sobrepico (fenómeno de Gibbs):**
- **El escalón perfecto solo se logra con la sinc infinita.** Nosotros la cortamos de golpe a 101 valores, que es multiplicarla por una ventana rectangular.
- **Ese corte brusco en el tiempo produce oscilaciones en frecuencia,** que se amontonan donde el escalón es discontinuo (en el corte) y dejan un sobrepico de ~10 %.
- **Agrandar Nh no lo arregla.** Medido con Nh = 101, 201 y 401, el sobrepico da 10.0 %, 9.5 % y 9.2 %: las oscilaciones se angostan, pero la altura queda en ~9 %.

**Fila 3, izquierda, "y[n] = IFFT(X[k] H[k])":** la misma salida que en la Figure 1 de `ej2` (cosenoidal limpia de amplitud 10 con transitorios, de n = −50 a 450) y después **todo ceros**. Que el final dé 0 es la prueba de que N = 1024 alcanzó: si N fuera menor que 501, el final de y se "doblaría" sobre el principio (convolución circular) y lo arruinaría.

**Fila 3, derecha, "|Y[k]| = |X[k]| |H[k]|":**
- **Se mantiene el pico de ≈ 1922 en k = 15** (y su espejo en 1009).
- **El pico de k = 102 desapareció.**

**Por qué:** en cada k, |Y| = |X|·|H|. En k = 15, |H| ≈ 1, así que |Y| ≈ |X| y el pico queda igual. En k = 102, |H| ≈ 0, así que |Y| ≈ 0 y el pico se elimina. **Es el mismo filtrado del Ej. 2, visto en frecuencia.**

### Conclusiones

- La convolución por FFT da **el mismo resultado** que `conv`. Para señales largas es mucho más rápida: ~N·log N operaciones contra ~Nx·Nh.
- Hace falta el **relleno de ceros** para evitar la convolución circular.
- **El sobrepico de Gibbs en |H| viene de truncar la sinc** y no desaparece agrandando Nh. Por eso el Ej. 4 usa otra ventana.
- **Entrada vs salida:** comparando |X| (fila 1) con |Y| (fila 3) se ve qué frecuencias pasaron (150 Hz) y cuáles se eliminaron (1000 Hz).

---

## 6. Ejercicio 4: filtro FIR con ventana de Hamming

**Archivos:** [`ej4.m`](ej4.m) (función) y [`ej4_estudios.m`](ej4_estudios.m)

| Se corre con | Abre | Qué muestra | PNG |
|---|---|---|---|
| `close all; ej4` | Figure 1 | El filtro con Hamming, 4×2 gráficos | `ej4_1.png` |
| `ej4_estudios` | Figure 1 | a) Rectangular vs Hamming | `ej4_estudios_1.png` |
| | Figure 2 | b) Banda de transición vs Nh | `ej4_estudios_2.png` |
| | Figure 3 | c) Espectro de entrada vs Nx | `ej4_estudios_3.png` |
| | Figure 4 | d) Variación de Fm | `ej4_estudios_4.png` |

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

Valores: Fm = 10000 Hz, Nh = 101, Fc = 500 Hz, Nx = 401, N = 1024.

**Por qué normalizar:** la ventana cambia la suma de h (acá pasa a 1.0024), y esa suma es la ganancia en continua (f = 0). Dividiendo por ella, las bajas frecuencias pasan con amplitud exacta.

**Comparación conv vs FFT:** la consola imprime `max|y_fft - y_conv| = 1.066e-14`. Son idénticas.

---

### Figure 1 de `ej4`: "Ejercicio 4 - FIR Hamming (Fm = 10000 Hz, Nh = 101)"

Izquierda: tiempo en segundos, de **−0.005 a 0.045 s**. Derecha: frecuencia en Hz, de **0 a 5000 Hz** (Nyquist; no se muestra el espejo).

**Fila 1, izquierda, "Entrada x(t_n) (F_m = 10000 Hz, N_x = 401)":** la misma entrada de siempre, ahora con el eje en segundos. Va de 0 a 0.04 s (401 × 0.0001 s), y las ondas grandes de 150 Hz duran 0.00667 s cada una (se ven 6).

**Fila 1, derecha, "|X(f)|":** picos de ≈ 1922 en **150 Hz** y ≈ 962 en **1000 Hz**. Son los mismos de la Figure 1 de `ej3`, pero el eje ahora dice directamente los Hz (cada k se multiplicó por Fm/N = 9.77 Hz).

**Fila 2, izquierda, "h(t_n): sinc x Hamming normalizada (N_h = 101, F_c = 500 Hz)":**
- La sinc entre −0.005 y 0.005 s (n = −50 a 50), con máximo ≈ 0.1 en t = 0.
- Los lóbulos laterales (≈ −0.02) **se apagan suavemente hasta 0 en los extremos**.
- Comparala con el gráfico del medio de la Figure 1 de `ej2`: ahí la sinc terminaba de golpe en n = ±50 con lóbulos de ±0.007. Esa es la diferencia que introduce la ventana de Hamming.

**Fila 2, derecha, "|H(f)|":**
- **Plana en 1 hasta ≈ 350 Hz,** baja suavemente y **llega a 0 cerca de 650 Hz**.
- **Sin ondulaciones ni sobrepico:** el máximo es 1.0. En la Figure 1 de `ej3` el máximo era ≈ 1.1.

**Por qué no hay Gibbs:** como la sinc ya no se corta de golpe sino que se apaga gradualmente, desaparecen las oscilaciones en frecuencia. A cambio, la caída de 1 a 0 es menos abrupta: la **banda de transición** (de ≈ 350 a ≈ 650 Hz) es más ancha que con la sinc cortada.

**Filas 3 y 4, izquierda, "Salida y(t_n) por convolución en el tiempo (conv)" y "... por convolución FFT (ifft)":**
- Las dos son **idénticas**: una cosenoidal limpia de 150 Hz y amplitud 10, de 0 a 0.04 s, con transitorios de −0.005 a 0.005 s y de 0.035 a 0.045 s.
- Comparadas con el y del Ej. 2, los transitorios casi no tienen ondulaciones.

**Por qué son iguales:** la convolución en el tiempo y la multiplicación en frecuencia son matemáticamente equivalentes (sección 2.4). La diferencia numérica es de 10⁻¹⁴, solo redondeo.

**Filas 3 y 4, derecha, "|Y(f)| = |X(f)| |H(f)| (FFT)":** el **mismo gráfico** dos veces, como pide la consigna. Queda solo el pico de ≈ 1922 en 150 Hz; el de 1000 Hz desapareció, porque ahí |H| = 0.

---

### Figure 1 de `ej4_estudios`: "a) Rectangular vs Hamming"

Consigna: *"Verificar la mejora en el espectro H(f) reduciendo el fenómeno de Gibbs"*. Compara los dos núcleos con el mismo Nh = 101:
- **Curva azul, "Rectangular (Ej. 2 y 3)":** la sinc cortada de golpe.
- **Curva naranja, "Hamming normalizada (Ej. 4)".**

**Arriba, "|H(f)|, N_h = 101, F_c = 500 Hz" (escala lineal, de 0 a 5000 Hz):**
- **Azul:** ondula entre ≈ 0.95 y ≈ 1.05 de 0 a 350 Hz, tiene un **sobrepico de ≈ 1.1 cerca de 400 Hz**, cae a 0 cerca de 560 Hz y después muestra **pequeños rebotes** (lóbulos) hasta ≈ 1000 Hz.
- **Naranja:** **plana en 1**, sin sobrepico. Cae más suave (de ≈ 350 a ≈ 650 Hz) y después queda en 0, sin rebotes visibles.
- **Las dos cruzan 0.5 en 500 Hz,** que es la Fc.

**Abajo, "|H(f)| en dB":** permite ver los rebotes chicos, que en lineal parecen 0 (ver la tabla de dB en la sección 2.6).
- **Azul:** el primer lóbulo después del corte llega a ≈ **−22 dB** (≈ 0.08) cerca de 600 Hz. Los siguientes bajan despacio y siguen en ≈ −60 dB (≈ 0.001) a 5000 Hz.
- **Naranja:** después del corte, los lóbulos quedan entre ≈ −50 y ≈ −85 dB, mucho más abajo en toda la banda.

**Números que imprime la consola** (lóbulo máximo medido de 650 Hz para arriba):

| Núcleo (Nh = 101) | Sobrepico | Lóbulo lateral máx. (f > 650 Hz) |
|---|---|---|
| Rectangular | 10.0 % | −27.9 dB |
| Hamming | **0.0 %** | **−41.0 dB** |

(El −41 dB de Hamming se mide justo en 650 Hz, donde todavía está terminando de caer. Más lejos del corte, sus lóbulos quedan en −50 a −85 dB, como se ve en el gráfico.)

**Por qué:** la ventana de Hamming apaga los extremos de la sinc en vez de cortarlos de golpe. Sin corte brusco no hay oscilaciones fuertes en frecuencia: desaparece el sobrepico y los lóbulos bajan.

**Lo que se pierde:** la caída de la naranja es más lenta que la de la azul. Hamming tiene una **banda de transición más ancha**.

**Conclusión: Gibbs reducido, verificado.**

---

### Figure 2 de `ej4_estudios`: "b) Banda de transicion vs Nh"

Consigna: *"Investigar los cambios en la calidad de su comportamiento en la banda de transición ... al modificar las dimensiones del núcleo del filtro"*.

Un solo gráfico, **"|H(f)| con ventana de Hamming para distintos N_h"**, de 0 a 2000 Hz. Colores según la leyenda: **N_h = 21** azul, **51** naranja, **101** amarillo, **201** violeta.

**Qué se ve:**
- **Todas valen 1 cerca de 0 Hz** (por la normalización) y bajan alrededor de 500 Hz.
- **Azul (Nh = 21):** baja muy despacio. Empieza a caer desde ≈ 100 Hz, vale ≈ 0.6 en 500 Hz y recién llega a ≈ 0 cerca de 1350 Hz. **En 1000 Hz todavía vale ≈ 0.09.**
- **Naranja (51) y amarillo (101):** cada vez más empinadas.
- **Violeta (Nh = 201):** casi vertical. Plana en 1 hasta ≈ 450 Hz y en 0 desde ≈ 550 Hz.

**Números que imprime la consola** (transición medida entre |H| = 0.9 y |H| = 0.1):

| Nh | Transición | Ancho | \|H(1000 Hz)\| |
|---|---|---|---|
| 21 | 228 – 989 Hz | 760 Hz | 0.0933 |
| 51 | 316 – 681 Hz | 365 Hz | 0.0017 |
| 101 | 409 – 592 Hz | 183 Hz | 0.0007 |
| 201 | 455 – 546 Hz | 90 Hz | 0.0005 |

**Por qué:** un núcleo más largo incluye más lóbulos de la sinc ideal y se parece más a ella, así que su |H| se parece más al escalón ideal. El ancho de la transición es **inversamente proporcional a Nh**: duplicar Nh lo reduce a la mitad (Nh × ancho ≈ 18000).

**Consecuencia práctica:** con Nh = 21 la componente de 1000 Hz pasaría con 0.0933 × 5 ≈ 0.47 de amplitud. Con Nh ≥ 51 queda en menos de 0.01. El costo de un Nh grande es más cálculo y transitorios más largos.

---

### Figure 3 de `ej4_estudios`: "c) Efecto de enventanado sobre X(f)"

Consigna: *"Investigar los cambios en la calidad del espectro de la señal de entrada, al modificar el tamaño de la ventana no nula, y verificar el efecto enventanado comparándolo con el espectro de la señal de longitud infinita"*.

Tres gráficos, uno debajo del otro, de 0 a 1500 Hz. En cada uno:
- **Los palitos rojos son el espectro de la señal infinita:** dos rayas perfectas, de altura **10 en 150 Hz** y **5 en 1000 Hz** (las amplitudes del coseno y del seno).
- **La curva azul es 2|X(f)|/N_x de la señal cortada a N_x muestras.** Se divide por Nx y se multiplica por 2 para que la altura de cada pico sea comparable con la amplitud.

**Arriba, "2|X(f)|/N_x con N_x = 101 (ventana de 10.1 ms)":**
- Dos **lóbulos muy anchos**, de ≈ 200 Hz de base. El de 150 Hz llega a ≈ 10.5 (corrido un poco a ≈ 160 Hz) y el de 1000 Hz a ≈ 5.6.
- Entre los dos picos hay **muchos lóbulos laterales** de altura 1 a 3, separados ≈ 99 Hz (10000/101).

**Medio, "... N_x = 401 (ventana de 40.1 ms)":** los lóbulos principales se angostan a ≈ 50 Hz de base. Los picos coinciden con los rojos (≈ 10 y ≈ 5), y los lóbulos laterales son más chicos y más juntos.

**Abajo, "... N_x = 1601 (ventana de 160.1 ms)":** los picos son **casi líneas**, de ≈ 12 Hz de base, y prácticamente coinciden con las rayas rojas.

**Por qué:**
- **Una senoidal infinita tiene una sola frecuencia exacta,** que da una raya.
- **Tomar solo Nx muestras es multiplicarla por una ventana rectangular,** y en frecuencia cada raya se "desparrama" en un lóbulo de ancho ≈ 2·Fm/Nx: 198 Hz, 50 Hz y 12.5 Hz para los tres casos.
- **Con Nx = 101 los lóbulos son tan anchos que se pisan entre sí** (y con los de las frecuencias negativas). Por eso las alturas no coinciden exactamente con 10 y 5.
- **Cuanto más larga la ventana, mejor se distinguen frecuencias cercanas:** mejor **resolución en frecuencia**.

**Detalle extra:** al alargar la ventana, los lóbulos laterales se angostan pero su altura **relativa** al pico no baja. Siguen siendo ≈ 22 % del pico, porque es la misma ventana rectangular. Es el mismo fenómeno de Gibbs, ahora aplicado a la señal de entrada. Para bajarlos habría que aplicar una ventana suave (por ejemplo, Hamming) también a x.

---

### Figure 4 de `ej4_estudios`: "d) Variacion de Fm"

Consigna: *"Verificar las modificaciones en la calidad de representación de las señales en tiempo discreto y frecuencias al modificar la frecuencia de muestreo Fm"*.

Cuatro filas, una por cada Fm. En todas se mantiene fija, en segundos, la duración de la ventana (0.04 s) y del núcleo (≈ 0.01 s); solo cambia cuántas muestras se toman.

- **Columna izquierda:** x en azul y la salida y en naranja (el título dice "rojo", pero Octave lo dibuja naranja-rojizo), de 0 a 0.04 s.
- **Columna derecha:** |X(f)| de 0 a Nyquist.

| Fila | Fm | Nyquist | Nx | Nh |
|---|---|---|---|---|
| 1 | 10000 Hz | 5000 Hz | 401 | 101 |
| 2 | 3000 Hz | 1500 Hz | 121 | 31 |
| 3 | 2000 Hz | 1000 Hz | 81 | 21 |
| 4 | 1200 Hz | 600 Hz | 49 | 13 |

**Fila 1, "F_m = 10000 Hz":**
- **Izquierda:** x azul con la ondulación de 1000 Hz bien dibujada (10 muestras por ciclo). La naranja es la cosenoidal de 150 Hz suave, de amplitud 10. **El filtro funciona.**
- **Derecha:** picos en 150 Hz (≈ 1922) y en 1000 Hz (≈ 962).

**Fila 2, "F_m = 3000 Hz":**
- **Izquierda:** x azul se ve **quebrada, en zigzag**: la componente de 1000 Hz tiene solo 3 muestras por ciclo. La naranja sigue siendo la cosenoidal de 150 Hz limpia.
- **Derecha, de 0 a 1500 Hz:** los dos picos siguen en su lugar, 150 Hz (≈ 600) y 1000 Hz (≈ 300). Son más bajos porque hay menos muestras: 10·121/2 = 605 y 5·121/2 = 302.
- **Por qué funciona:** 1000 Hz < Nyquist = 1500 Hz. La componente está bien representada, aunque con pocas muestras se vea fea en el tiempo.

**Fila 3, "F_m = 2000 Hz":**
- **Izquierda:** **la azul y la naranja se superponen.** x ya es una cosenoidal pura de 150 Hz y amplitud 10, sin ondulación, y el filtro la deja igual.
- **Derecha, de 0 a 1000 Hz:** **un solo pico**, en 150 Hz (≈ 400 = 10·81/2). No hay nada en 1000 Hz.
- **Por qué:** 1000 Hz es justo Nyquist. El seno se muestrea en sus cruces por cero, 5·sen(π·n) = 0, y desaparece (lo mismo que en la Figure 1 de `ej2_variantes`).

**Fila 4, "F_m = 1200 Hz":**
- **Izquierda:** la azul y la naranja **se superponen otra vez**, pero ahora la forma es **irregular**: picos de distinta altura, de hasta ≈ 14.
- **Derecha, de 0 a 600 Hz:** **dos picos**, uno en 150 Hz (≈ 245) y **otro en 200 Hz** (≈ 120).
- **Por qué aparece un pico en 200 Hz:** es **aliasing**. 1000 Hz > Nyquist = 600 Hz, así que la componente aparece disfrazada en 1200 − 1000 = **200 Hz**. Su pico mide la mitad que el de 150 Hz porque su amplitud es 5 (≈ 5·49/2 = 122).
- **Por qué y = x (el filtro no la elimina):** el alias está en 200 Hz, **por debajo de Fc = 500 Hz**. Para el filtro es una frecuencia que tiene que pasar. La salida queda **contaminada** con una componente que en la señal original era de 1000 Hz, y la forma irregular es la suma de 150 Hz y 200 Hz.

**Conclusión:** **Fm tiene que superar el doble de la máxima frecuencia de la señal**, acá 2 × 1000 = 2000 Hz.
- **Con Fm de sobra (fila 1)** todo se representa bien.
- **Con Fm justa (fila 2)** se representa, pero con pocas muestras.
- **En el límite (fila 3)** se puede perder la componente.
- **Por debajo (fila 4)** aparece como una frecuencia falsa que ningún filtro digital puede separar. La única solución es filtrar en analógico *antes* de muestrear (filtro **antialiasing**).

---

## 7. Diferencias con el enunciado

| Enunciado | En el TP | Motivo |
|---|---|---|
| "frecuencias de 30 Hz y 1000 Hz" | 150 Hz | La fórmula dice cos(2π·**150**·t). Error de tipeo de la guía |
| Ejemplo Fc = 1000 Hz | Fc = 500 Hz (1000 Hz se explica aparte) | La consigna pide Fc *entre* las dos frecuencias. Con 1000 Hz, esa componente pasa a la mitad |
| Comando `haming` | `hamming` | Error de tipeo de la guía |
| "su correspondiente normalizado" | `h / sum(h)` | Interpretación habitual: ganancia unitaria en continua |
| Espectros en Hz | Se muestran de 0 a Fm/2 | La otra mitad es espejo (simetría de la FFT de señales reales) |

---

## 8. Cómo lo toman: guion de defensa ejercicio por ejercicio

Según el relato de otro grupo, el docente va ejercicio por ejercicio:

| Ej. | Qué pidió |
|---|---|
| 1 | Explicar cómo se hace la convolución y cómo se llega al resultado |
| 2 | Identificar el cos y el sen. Después hizo **cambiar Δt y Fc en vivo, por separado**, y preguntó qué señal resultaba |
| 3 | Lo mismo pero en frecuencia: reconocer el **fenómeno de Gibbs** |
| 4 | Que es una mejora del 3: con la **ventana de Hamming** la señal queda más limpia |

Para cada uno, abajo hay qué decir, cómo hacerlo en Octave y las trampas posibles. Las tablas del Ej. 2 salen de simular cada caso en Octave 11.3.

### 8.1 Ejercicio 1: explicar la convolución

Con la Figure 1 de `ej1` abierta.

**Qué decir (en 30 segundos):**

> "La salida de un sistema LTI es la convolución de la entrada con la respuesta al impulso: y[n] = Σₖ x[k]·h[n−k]. Para cada n, se da vuelta h, se la corre n lugares, se multiplica con x donde se superponen y se suma. Como x es un pulso de unos, cada y[n] es la suma de los valores de la rampa que caen dentro del pulso. Primero crece mientras la rampa entra, llega al máximo 15 cuando entra entera, y después baja porque la rampa termina y se van perdiendo valores. El largo de y es 6 + 6 − 1 = 11 y arranca en n = −3 + 0."

La tabla de la ventana deslizante está en la sección 3.

**Trampas y cómo responderlas:**

- **"¿Por qué hay salida donde x vale 0 (n = 3 a 7)?"** Porque y[n] no es x[n]·h[n]. Depende de valores de x **anteriores** a n: el sistema tiene "memoria", como el eco en una sala.
- **"¿Por qué y vale 0 en −3 aunque x[−3] = 1?"** Porque h[0] = 0, así que y[−3] = x[−3]·h[0] = 0.
- **"¿Por qué baja si x es todo unos y h sube?"** Porque las dos son **finitas**, como dice el enunciado ("en un intervalo de tiempo guardado en el vector nx/nh"). Cuando la rampa termina, ya no entra nada nuevo. Si fueran infinitas, y crecería sin límite: (n+3)(n+4)/2.
- **"¿Cuándo la salida es igual a h?"** Solo si x = δ[n] (un único 1 en n = 0). Un pulso es una suma de impulsos corridos, así que la salida es una suma de copias corridas de h.

### 8.2 Ejercicio 2: identificar la señal y cambiar Δt y Fc

#### Identificar el cos y el sen

Con la Figure 1 de `ej2` abierta. En x(t) = **10** cos(2π·**150**·t) + **5** sen(2π·**1000**·t):

- **Amplitud** = el número de adelante (10 y 5).
- **Frecuencia** = el número que multiplica a 2π·t (150 Hz y 1000 Hz).
- **En el gráfico de x[n]:** la onda lenta y grande es el coseno de 150 Hz (6 ciclos) y la ondulación rápida montada encima es el seno de 1000 Hz (40 ciclos, 10 muestras cada uno).
- **x[0] = 10:** el coseno arranca en su máximo y el seno en 0.

**Qué hace el sistema:** es un pasa bajos (sinc truncada) con corte en Fc. Con Fc = 500 Hz deja pasar 150 Hz y elimina 1000 Hz, así que **la señal resultante es ≈ 10 cos(2π·150·t)**, con transitorios en los bordes.

#### Cómo cambiar los parámetros en vivo

La función es `ej2(dt, Nh, Fc)`. Un `[]` deja el valor por defecto. Conviene poner `close all;` antes para que la ventana nueva sea la Figure 1:

```octave
close all; ej2(0.0005)            % cambia solo dt
close all; ej2([], [], 300)       % cambia solo Fc (dt = 0.0001, Nh = 101)
close all; ej2(0.0001, 101, 1000) % todo explícito
close all; ej2([], 21)            % cambia solo Nh
```

#### Procedimiento para predecir la salida (hacerlo en voz alta)

1. **Calcular Fm = 1/Δt y Nyquist = Fm/2.**
2. **Revisar cada componente de x.**
   - Si su frecuencia es menor que Nyquist, se representa bien.
   - Si es mayor, hay **aliasing** y aparece en |f − k·Fm| (el más cercano).
   - Si es un **seno** justo en Nyquist o en un múltiplo de Fm, sus muestras valen 0 y **desaparece**.
3. **Comparar la frecuencia (ya con alias, si lo hubo) con Fc.**
   - Debajo de Fc pasa con amplitud ≈ igual.
   - Encima de Fc se elimina.
   - Cerca de Fc pasa **a la mitad**.
4. **Revisar 2·Fc·Δt.** Si da 1, h[n] = sinc(n) = δ[n] y el filtro no hace nada (y = x). Si da más de 1, Fc quedó por encima de Nyquist y la fórmula deja de ser un pasa bajos.

#### Tabla: variando Fc (Δt = 0.0001 s, Nh = 101)

Amplitudes medidas en la zona central de y (sin transitorios):

| Fc | Amp. 150 Hz | Amp. 1000 Hz | Señal resultante |
|---|---|---|---|
| 50 Hz | 1.4 | 0.0 | Casi nada: los 150 Hz se eliminan (queda un resto porque con Fc tan baja la sinc es muy ancha y Nh = 101 la corta mucho) |
| 100 Hz | 0.6 | 0.0 | Casi nula |
| **150 Hz** | **5.3** | 0.0 | ≈ 5 cos(2π·150t): Fc **justo en** 150 Hz, pasa **a la mitad** |
| 300 Hz | 10.1 | 0.0 | ≈ 10 cos(2π·150t) |
| **500 Hz** (TP) | 10.0 | 0.07 | ≈ 10 cos(2π·150t) |
| 800 Hz | 10.0 | 0.2 | ≈ 10 cos(2π·150t), con un resto mínimo de 1000 Hz |
| **1000 Hz** | 10.0 | **2.5** | 10 cos(2π·150t) + **2.5** sen(2π·1000t): el seno pasa **a la mitad** |
| 1200 Hz | 10.0 | 4.7 | ≈ x(t): pasan casi enteras las dos |
| 2000 Hz o más | 10.0 | 5.0 | = x(t): el filtro no elimina nada |

**Frase para cerrar:** *"Lo que está por debajo de Fc pasa, lo que está por encima se elimina, y lo que cae justo en Fc pasa a la mitad."*

#### Tabla: variando Δt (Fc = 500 Hz, Nh = 101)

| Δt | Fm | Nyquist | 2FcΔt | Amp. 150 | Amp. 1000 | Señal resultante y por qué |
|---|---|---|---|---|---|---|
| 0.00005 s | 20000 Hz | 10000 Hz | 0.05 | 9.9 | 0 | 10 cos(2π·150t). Filtra bien; las 401 muestras duran solo 0.02 s |
| **0.0001 s** | 10000 Hz | 5000 Hz | 0.1 | 10.0 | 0.07 | 10 cos(2π·150t). Caso del TP |
| 0.0002 s | 5000 Hz | 2500 Hz | 0.2 | 10.2 | 0.04 | 10 cos(2π·150t) |
| 0.0003 s | 3333 Hz | 1667 Hz | 0.3 | 10.0 | 0.03 | 10 cos(2π·150t). Los 1000 Hz tienen solo 3.3 muestras por ciclo, pero se representan y se filtran |
| **0.0005 s** | 2000 Hz | **1000 Hz** | 0.5 | 10.0 | **0** | 10 cos(2π·150t), pero porque el seno de 1000 Hz cae **en Nyquist**: sen(π·n) = 0, desaparece **al muestrear** |
| **0.0008 s** | 1250 Hz | 625 Hz | 0.8 | 10.0 | **5.0** | **Aliasing:** los 1000 Hz aparecen en 1250 − 1000 = **250 Hz**. Como 250 < Fc, **pasa entero**: salida = 150 Hz + alias de 250 Hz |
| **0.001 s** | 1000 Hz | 500 Hz | **1.0** | 10.0 | 0 | 10 cos(2π·150t). Aliasing: sen(2π·n) = 0. Además h = δ, así que **y = x**: el filtro no hace nada |
| **0.002 s** | 500 Hz | 250 Hz | **2.0** | **20.0** | 0 | **20** cos(2π·150t). h[n] = 2·sinc(2n) = 2δ[n], o sea **y = 2x**. Fc = 500 Hz quedó por encima de Nyquist y la fórmula ya no es un pasa bajos |

**Frase para cerrar:** *"Si Δt es chico (Fm alta) todo se representa bien y el filtro trabaja. Si Δt crece y Fm baja de 2000 Hz, los 1000 Hz no cumplen Nyquist: o desaparecen o aparecen como otra frecuencia (alias), y si el alias queda debajo de Fc, el filtro no lo puede sacar."*

**Por qué cambia el gráfico aunque la señal resultante sea la misma:** con las mismas 401 muestras, un Δt más grande abarca más tiempo (Nx·Δt) y se ven más ciclos, con menos muestras por ciclo. Compará la Figure 1 de `ej2` con las Figures 1 y 2 de `ej2_variantes`.

### 8.3 Ejercicio 3: frecuencia y fenómeno de Gibbs

Con la Figure 1 de `ej3` abierta.

**Qué decir:**

> "Es el mismo filtrado del Ej. 2, pero hecho en frecuencia: la convolución en el tiempo equivale a multiplicar las transformadas. Se completan x y h con ceros hasta N = 1024 para que la convolución circular de la FFT coincida con la lineal (hace falta N ≥ 401 + 101 − 1 = 501). Se calculan las FFT, se multiplican y se antitransforma. El resultado es igual al de conv (diferencia de 10⁻¹⁴, solo redondeo)."

**Cómo leer los espectros (k → f = k·Fm/N ≈ k·9.77 Hz):**

- **|X[k]| (fila 1):** dos picos, en k ≈ 15 (150 Hz) y en k ≈ 102 (1000 Hz), más sus espejos del otro lado de k = 512.
- **|H[k]| (fila 2):** un escalón. Vale ≈ 1 hasta k ≈ 51 (Fc = 500 Hz) y ≈ 0 después, más su espejo.
- **|Y[k]| = |X|·|H| (fila 3):** queda solo el pico de 150 Hz.

**Fenómeno de Gibbs, qué decir exactamente:**

> "En |H[k]| aparecen ondulaciones y un sobrepico de ~10 % justo antes del salto de 1 a 0, en el corte. Se debe a que la sinc ideal es infinita y la cortamos de golpe a 101 valores, que es como multiplicarla por una ventana rectangular. Ese corte abrupto en el tiempo genera las oscilaciones en frecuencia, que se concentran donde el filtro ideal es discontinuo. Aunque se agrande Nh no desaparece: el sobrepico se mantiene en ~9 %, solo se angosta."

Ojo con la frase "por el corte abrupto en el medio del gráfico". La **causa** es el corte abrupto de la sinc en el **tiempo**. Las oscilaciones **aparecen** en el salto de |H| (en Fc), que en el gráfico de 0 a 1023 se ve dos veces: cerca de k ≈ 51 y de su espejo, k ≈ 973.

**Si te hacen cambiar Fc acá** (`close all; ej3([], [], 1000)`): el escalón de |H[k]| se mueve a k ≈ Fc·N·Δt. Con Fc = 1000 Hz, k ≈ 102, justo sobre el pico de 1000 Hz, que queda a la mitad. Gibbs sigue estando en el nuevo corte.

### 8.4 Ejercicio 4: Hamming

Con la Figure 1 de `ej4` abierta (y, si hay tiempo, la Figure 1 de `ej4_estudios`).

**Qué decir:**

> "Es una mejora del Ej. 3. Se multiplica la sinc por una ventana de Hamming, una campana que apaga suavemente los extremos en vez de cortarlos de golpe. Así se reduce el fenómeno de Gibbs: |H(f)| ya no tiene sobrepico (máximo 1.0 contra 1.1) y los lóbulos laterales bajan (de −28 dB a −41 dB o menos). Después se normaliza dividiendo por la suma de h, para que la ganancia en bajas frecuencias sea 1. El precio es que la transición entre la banda que pasa y la eliminada es un poco más ancha. Además se usa Fm como dato en vez de Δt, y los ejes quedan en segundos y en Hz. La salida por conv y por FFT es la misma."

**Si preguntan algo más:**

- **Nh más grande** (Figure 2 de `ej4_estudios`): transición más angosta (se reduce a la mitad al duplicar Nh).
- **Nx más grande** (Figure 3): picos de |X| más finos (mejor resolución), más parecidos a las rayas de la señal infinita.
- **Fm más baja** (Figure 4): con menos de 2000 Hz aparece aliasing de la componente de 1000 Hz.

### 8.5 Antes de la defensa

- [ ] Tener Octave abierto **en la carpeta del TP** y probar `ej1`, `ej2`, `ej3` y `ej4`.
- [ ] Practicar `close all; ej2([], [], 1000)` y `close all; ej2(0.0008)` y saber explicar qué se ve.
- [ ] Saber hacer a mano y[−1] o y[0] del Ej. 1.
- [ ] Saber pasar de k a Hz: f = k·Fm/N.
- [ ] Tener a mano la regla de Nyquist: Fm > 2·f_max.

---

## 9. Preguntas para preparar la defensa

**¿Por qué hay que calcular `ny` a mano?**
Porque `conv` solo devuelve valores; no sabe en qué n empieza cada señal. Inicio de y = inicio de x + inicio de h.

**¿Por qué Nh tiene que ser impar?**
Para tener una muestra en n = 0 y la misma cantidad a cada lado: h queda simétrica y centrada, y la salida no se retrasa.

**¿Qué significa 2·Fc·Δt en la fórmula de h?**
Es la frecuencia de corte normalizada (Fc dividida por Nyquist). Define el valor de h[0] y cada cuántas muestras cruza por cero la sinc (cada 1/(2FcΔt)).

**¿Por qué con Δt = 0.001 s el filtro "no hace nada"?**
Porque 2·Fc·Δt = 1 y sinc(n) es 0 en todo entero salvo n = 0, así que h = δ e y = x. Además Fc = 500 Hz coincide con Nyquist.

**¿Por qué desaparece la componente de 1000 Hz con Δt = 0.0005 s?**
Porque 1000 Hz es justo Fm/2. Las muestras de un seno en Nyquist caen todas en sus cruces por cero.

**¿Por qué con Nh = 201 la amplitud da 10.2 y no 10?**
Porque la sinc sigue cortada de golpe y la ganancia ondula ±5 % en la banda de paso (Gibbs). Con Nh = 201, los 150 Hz caen en una cresta de esa ondulación.

**¿Por qué hay que rellenar con ceros antes de la FFT?**
La DFT supone señales periódicas, y el producto X·H da una convolución circular. Con N ≥ Nx + Nh − 1 no hay solapamiento y coincide con la lineal.

**¿Por qué N potencia de 2?**
Porque así el algoritmo FFT es más eficiente.

**¿Por qué el pico de 150 Hz en |X[k]| mide el doble que el de 1000 Hz?**
Porque la altura es ≈ amplitud·Nx/2, y la amplitud del coseno (10) es el doble que la del seno (5).

**¿Qué es el fenómeno de Gibbs y de dónde sale?**
Son las ondulaciones y el sobrepico (~9 %) que aparecen en |H| cerca del corte por truncar la sinc con una ventana rectangular. Agrandar Nh no lo elimina.

**¿Qué mejora la ventana de Hamming y qué empeora?**
Elimina el sobrepico y baja los lóbulos laterales (de −28 dB a −41 dB o menos). A cambio, la transición es más ancha.

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
Una frecuencia mayor a Fm/2 que aparece como otra más baja. Por ejemplo, 1000 Hz con Fm = 1200 Hz aparece en 200 Hz y pasa el filtro (Figure 4 de `ej4_estudios`, fila 4).

**¿Cuál es la diferencia entre `.*` y `*` en Octave?**
`.*` multiplica elemento a elemento; `*` es el producto matricial.

**¿Por qué `real(ifft(...))`?**
Porque el redondeo deja una parte imaginaria de ~10⁻¹⁵ que no tiene significado físico.

---

## 10. Mapa de archivos y figuras

| Archivo | Tipo | Contenido |
|---|---|---|
| [`ej1.m`](ej1.m) | script | Ejercicio 1 |
| [`ej2.m`](ej2.m) | función `ej2(dt, Nh, Fc, Nx)` | Ejercicio 2 (un caso) |
| [`ej2_variantes.m`](ej2_variantes.m) | script | Las variantes de Δt y Nh (4 gráficos) |
| [`ej3.m`](ej3.m) | función `ej3(dt, Nh, Fc, Nx, N)` | Ejercicio 3 |
| [`ej4.m`](ej4.m) | función `ej4(Fm, Nh, Fc, Nx, N)` | Ejercicio 4 |
| [`ej4_estudios.m`](ej4_estudios.m) | script | Estudios a) a d) del Ej. 4 |
| [`correr_todo.m`](correr_todo.m) | script | Corre todo y guarda los PNG |
| [`guardar_figuras.m`](guardar_figuras.m) | función auxiliar | Guarda las figuras abiertas |
| `figuras_octave/` | imágenes | Figuras generadas |

**Qué figura es cada PNG:**

| PNG | Ventana de Octave | Contenido |
|---|---|---|
| `ej1_1.png` | Figure 1 de `ej1` | Convolución pulso × rampa |
| `ej2_1.png` | Figure 1 de `ej2` | Caso base: Δt = 0.0001 s, Nh = 101 |
| `ej2_variantes_1.png` | Figure 1 de `ej2_variantes` | Δt = 0.0005 s |
| `ej2_variantes_2.png` | Figure 2 de `ej2_variantes` | Δt = 0.001 s |
| `ej2_variantes_3.png` | Figure 3 de `ej2_variantes` | Nh = 21 |
| `ej2_variantes_4.png` | Figure 4 de `ej2_variantes` | Nh = 201 |
| `ej3_1.png` | Figure 1 de `ej3` | Convolución por FFT, 3×2 |
| `ej4_1.png` | Figure 1 de `ej4` | Filtro con Hamming, 4×2 |
| `ej4_estudios_1.png` | Figure 1 de `ej4_estudios` | a) Rectangular vs Hamming |
| `ej4_estudios_2.png` | Figure 2 de `ej4_estudios` | b) Transición vs Nh |
| `ej4_estudios_3.png` | Figure 3 de `ej4_estudios` | c) Espectro vs Nx |
| `ej4_estudios_4.png` | Figure 4 de `ej4_estudios` | d) Variación de Fm |

Todos los números de esta guía salen de correr los scripts en Octave 11.3.
