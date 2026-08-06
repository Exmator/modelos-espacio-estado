# Librerías de R para testing, validación de modelos y simulación

Este documento reúne las librerías de R relevantes para el proyecto: comparar
estadística y computacionalmente tres modelos de espacio-estado. Se agrupan en
tres categorías: testing unitario, validación de modelos de estimación, y
manejo de estudios de simulación.

---

## 1. Testing unitario

### `testthat`
El framework de testing estándar de facto en R (usado por la mayoría de
paquetes en CRAN). Permite escribir pruebas con la sintaxis
`test_that("descripción", { expect_equal(...) })`, con funciones como
`expect_equal(x, y, tolerance = 1e-8)`, `expect_error()`, `expect_true()`.

**Aplicación en el proyecto:** reemplazar los `cat()`/`print()` de
`algoritmo_estimacion_test.R` y `algoritmo_simulacion_test.R` por
aserciones automáticas. Por ejemplo, la función `probar_error()` que ya
existe en `algoritmo_simulacion_test.R` (línea 99) es esencialmente un
`expect_error()` hecho a mano — se puede convertir directamente. También
sirve para testear de forma exacta las partes deterministas del algoritmo
MMP (`compute_Omega_c()`, `mmp_precompute()`), que no dependen de números
aleatorios y por tanto no tienen el problema de "flakiness".

### `tinytest`
Alternativa a `testthat`, más liviana y sin dependencias externas.

**Aplicación en el proyecto:** opción si se prefiere no depender del
ecosistema `testthat`/`devtools`; funcionalmente cubre lo mismo para este
proyecto.

### `waldo`
Complementa a `testthat` con diffs legibles de objetos complejos (listas,
matrices, arrays anidados).

**Aplicación en el proyecto:** útil para comparar la lista completa que
devuelve `estimar_estados_mmp()` (`mu_mat`, `draws`) o `modelo_de_simulacion()`
(`Y`, `alpha`, `Z`, `T_mat`, ...) contra un resultado de referencia sin tener
que escribir una aserción por cada campo.

---

## 2. Validación de modelos de estimación / suavizado de estados

Estos paquetes no son "de testing" en sí, sino implementaciones maduras de
modelos de espacio-estado que sirven como **referencia externa** (ground
truth) contra la cual contrastar la implementación propia del método MMP.

### `KFAS`
Filtro y suavizador de Kalman completo, implementado en C++ (rápido), soporta
modelos de espacio-estado gaussianos generales, incluyendo modelos
factoriales dinámicos.

**Aplicación en el proyecto:** correr el mismo modelo factorial dinámico
(la aplicación concreta usada en McCausland et al. 2011) en `KFAS` y
verificar que `E[α|y]` coincide con `mu_mat` de `estimar_estados_mmp()`.
Es la validación cruzada más directa disponible.

### `FKF` (Fast Kalman Filter)
Implementación minimalista del filtro de Kalman, orientada a velocidad.

**Aplicación en el proyecto:** útil como segundo punto de comparación de
velocidad, ya que el artículo de McCausland compara explícitamente métodos
basados en el filtro de Kalman contra los métodos basados en precisión
(CFA, MMP) — comparación central de la tesis.

### `dlm`
Paquete clásico (Petris, Petrone & Campagnoli) para modelos lineales
dinámicos (*Dynamic Linear Models*), con `dlmSmooth()` para suavizamiento.

**Aplicación en el proyecto:** tercera alternativa de referencia para
contrastar resultados, con una API distinta a `KFAS` que puede servir para
detectar errores que sean específicos de un solo paquete de referencia.

### `MARSS`
*Multivariate Autoregressive State-Space* — pensado específicamente para
modelos factoriales/multivariados como el del proyecto (varias series
observadas explicadas por pocos factores latentes).

**Aplicación en el proyecto:** al estar diseñado para el mismo tipo de
modelo factorial dinámico que usa McCausland et al. (2011), puede servir
como un cuarto modelo de comparación, o como referencia adicional de
validación más cercana estructuralmente al caso de uso del proyecto.

---

## 3. Manejo de estudios de simulación

### `SimDesign`
Paquete diseñado específicamente para estudios Monte Carlo que comparan
métodos bajo distintas condiciones experimentales (*design factors*).
Organiza el flujo generar datos → estimar → evaluar → agregar resultados,
maneja errores de forma robusta y reporta métricas agregadas sobre muchas
repeticiones (sesgo, RMSE, cobertura de intervalos, tiempo de ejecución).

**Aplicación en el proyecto:** es el más relevante para el objetivo central
de la tesis. Los "Caso 1: m=4, p=10" y "Caso 2: m=10, p=100" de
`algoritmo_simulacion_test.R` son exactamente *design factors* en el
vocabulario de `SimDesign`. Permite automatizar la comparación de los tres
modelos de espacio-estado sobre múltiples repeticiones y condiciones,
reportando de forma sistemática tanto la calidad de la estimación de
estados como el costo computacional — en vez de correr casos sueltos
manualmente como se hace ahora.

### `simstudy`
Librería para definir y generar datos simulados de forma declarativa
(distribuciones, correlaciones, fórmulas).

**Aplicación en el proyecto:** menos central que `SimDesign` (que ya cubre
el ciclo completo de comparación), pero útil si se necesitan escenarios de
simulación más variados o complejos que los definidos manualmente en
`algoritmo_simulacion.R`.

### `microbenchmark` / `bench`
Miden tiempos de ejecución con alta precisión, ejecutando cada expresión
muchas veces y reportando estadísticas robustas (mediana, cuantiles) en vez
de una sola medición ruidosa.

**Aplicación en el proyecto:** eje "computacional" de la comparación de la
tesis — medir de forma rigurosa cuál de los tres modelos de espacio-estado
es más eficiente en tiempo de ejecución, replicando el tipo de comparación
que hace McCausland et al. (2011) entre el filtro de Kalman, CFA y MMP
(Sección 3.2 del artículo, "Computational experiments"). `bench::mark()`
además reporta uso de memoria, relevante para el caso "data-rich" (p grande).

### `future` / `furrr`
Permiten paralelizar cómputo en R de forma sencilla (`future_map()` como
versión paralela de `purrr::map()`).

**Aplicación en el proyecto:** si el estudio Monte Carlo crece (más
repeticiones, más combinaciones de `m`/`p`/`n`), paralelizar las
repeticiones reduce el tiempo total del estudio comparativo. `SimDesign`
tiene soporte nativo para paralelización usando estos paquetes por debajo.

---

## Nota: testing unitario con código estocástico

`modelo_de_simulacion()` y `mmp_draw()` usan `rnorm()`, por lo que un
`expect_equal` ingenuo puede fallar de forma intermitente aunque el código
esté correcto. Estrategias para evitarlo:

1. **`set.seed()` local** dentro de cada `test_that()` — vuelve determinista
   la secuencia aleatoria, útil para tests de regresión (¿el resultado sigue
   siendo el mismo que la última vez que se verificó como correcto?).
2. **Separar lo determinista de lo aleatorio** — `compute_Omega_c()` y
   `mmp_precompute()` no llaman `rnorm()` y se pueden testear con
   `expect_equal(..., tolerance = 1e-10)` sin ningún riesgo de flakiness.
3. **Inyectar el ruido en vez de generarlo internamente** — refactorizar
   `mmp_draw()` para recibir la matriz de innovaciones como parámetro
   opcional, de forma que el test le pase valores fijos conocidos.
4. **Tolerancias Monte Carlo amplias** (varias desviaciones estándar del
   error Monte Carlo) en vez del intervalo de confianza "de reporte" (ej.
   95%), que por diseño fallaría ~5% de las veces. `SimDesign` resuelve esto
   de forma nativa al agregar sobre muchas repeticiones y reportar tasas de
   cobertura en vez de un resultado binario ruidoso.

---

## Resumen de recomendación

| Necesidad | Librería recomendada |
|---|---|
| Tests unitarios generales | `testthat` |
| Validar corrección del suavizador MMP | `KFAS` (referencia externa) |
| Comparar los 3 modelos de espacio-estado (Monte Carlo) | `SimDesign` |
| Medir eficiencia computacional | `bench` o `microbenchmark` |
| Acelerar el estudio de simulación | `future` / `furrr` |
