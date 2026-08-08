# Fórmulas del modelo factorial dinámico
## Aplicación a McCausland, Miller y Pelletier (2011)

Este documento particulariza las fórmulas del apéndice de McCausland, Miller y Pelletier (2011) al modelo factorial dinámico lineal y gaussiano.

> **Nota de notación.** En el artículo, la letra `v_t` se utiliza para la innovación apilada. En el modelo factorial dinámico suele utilizarse `v_t` para el shock del estado. Para evitar confusiones, aquí se llamará `e_t` a la innovación apilada y `v_t` al shock del estado.

---

## 1. Modelo general de McCausland et al.

El modelo lineal general considerado en el artículo puede escribirse como

\[
y_t = X_t\beta + Z_t\alpha_t + G_t u_t,
\]

\[
\alpha_{t+1}=W_t\beta+T_t\alpha_t+H_tu_t,
\qquad t=1,\ldots,n-1,
\]

con distribución inicial

\[
\alpha_1\sim N(a_1,P_1),
\]

\[
u_t\sim N(0,I).
\]

Aquí:

- `y_t` es el vector observado.
- `alpha_t` es el vector de estados latentes.
- `Z_t` conecta los estados con las observaciones.
- `T_t` describe la evolución temporal de los estados.
- `G_t` y `H_t` distribuyen el ruido común `u_t` entre la ecuación de observación y la ecuación de estado.

McCausland et al. definen la innovación apilada como

\[
e_t=
\begin{pmatrix}
G_tu_t\\
H_tu_t
\end{pmatrix}.
\]

La matriz `A_t` es la matriz de precisión de `e_t`:

\[
A_t=\operatorname{Var}(e_t)^{-1}.
\]

---

## 2. Modelo factorial dinámico

El modelo factorial dinámico usado como aplicación en el artículo es

\[
y_t=Z\alpha_t+u_t,
\qquad u_t\sim N_p(0,D),
\]

\[
\alpha_1\sim N_m(a,Q_1),
\]

\[
\alpha_{t+1}=T\alpha_t+v_t,
\qquad v_t\sim N_m(0,Q),
\qquad t=1,\ldots,n-1.
\]

Las dimensiones son

\[
y_t\in\mathbb{R}^{p},
\qquad
\alpha_t\in\mathbb{R}^{m},
\]

\[
Z\in\mathbb{R}^{p\times m},
\qquad
D\in\mathbb{R}^{p\times p},
\]

\[
T,Q,Q_1\in\mathbb{R}^{m\times m},
\qquad
 a\in\mathbb{R}^{m}.
\]

`p` es el número de series observadas y `m` es el número de factores latentes.

---

## 3. Reducción del modelo general

Para obtener el modelo factorial dinámico a partir del modelo general de McCausland et al., se puede definir un ruido base

\[
\xi_t=
\begin{pmatrix}
 u_t\\
 v_t
\end{pmatrix},
\qquad
\xi_t\sim N\left(0,
\begin{pmatrix}
D&0\\
0&Q
\end{pmatrix}
\right).
\]

Las matrices de selección son

\[
G_t=\begin{pmatrix}I_p&0\end{pmatrix},
\qquad
H_t=\begin{pmatrix}0&I_m\end{pmatrix}.
\]

Por tanto,

\[
G_t\xi_t=u_t,
\qquad
H_t\xi_t=v_t.
\]

La innovación apilada es entonces

\[
e_t=
\begin{pmatrix}
G_t\xi_t\\
H_t\xi_t
\end{pmatrix}
=
\begin{pmatrix}
u_t\\v_t\end{pmatrix}.
\]

Para no confundir el error de observación `u_t` con el ruido base del artículo, en este documento se utiliza directamente la forma independiente anterior.

---

## 4. Matriz de precisión `A_t`

La covarianza de la innovación apilada es

\[
\operatorname{Var}(e_t)=
\begin{pmatrix}
D&0\\
0&Q
\end{pmatrix}.
\]

Como los errores de observación y de estado son independientes,

\[
A_t=\operatorname{Var}(e_t)^{-1}
=
\begin{pmatrix}
D^{-1}&0\\
0&Q^{-1}
\end{pmatrix}.
\]

Al particionar `A_t` en bloques,

\[
A_t=
\begin{pmatrix}
A_{11,t}&A_{12,t}\\
A_{21,t}&A_{22,t}
\end{pmatrix},
\]

se obtiene

\[
A_{11,t}=D^{-1},
\]

\[
A_{12,t}=A_{21,t}=0,
\]

\[
A_{22,t}=Q^{-1}.
\]

Si `D` y `Q` no cambian con el tiempo, entonces `A_t` tampoco cambia y puede escribirse simplemente como `A`.

---

## 5. Distribución condicional de los estados

Apilando todos los estados,

\[
\alpha=
\begin{pmatrix}
\alpha_1\\
\alpha_2\\
\vdots\\
\alpha_n
\end{pmatrix},
\qquad
\alpha\in\mathbb{R}^{nm},
\]

la distribución condicional de los estados dados los datos es

\[
\alpha\mid y\sim N\left(\Omega^{-1}c,\Omega^{-1}\right).
\]

La media condicional o estado suavizado es

\[
E(\alpha\mid y)=\Omega^{-1}c.
\]

Para simular una trayectoria completa de estados se puede generar

\[
\alpha=\Omega^{-1}c+L^{-\top}z,
\qquad
z\sim N(0,I),
\]

donde `L` es una factorización de Cholesky de `Omega`, por ejemplo

\[
\Omega=LL^\top.
\]

---

## 6. Estructura de `Omega`

La matriz de precisión tiene estructura tridiagonal por bloques:

\[
\Omega=
\begin{pmatrix}
\Omega_{11}&\Omega_{12}&0&\cdots&0\\
\Omega_{21}&\Omega_{22}&\Omega_{23}&\cdots&0\\
0&\Omega_{32}&\Omega_{33}&\cdots&0\\
\vdots&\vdots&\vdots&\ddots&\Omega_{n-1,n}\\
0&0&0&\Omega_{n,n-1}&\Omega_{nn}
\end{pmatrix}.
\]

Cada bloque `Omega_{ts}` tiene dimensión `m x m`.

La estructura tridiagonal aparece porque cada estado `alpha_t` está relacionado directamente solamente con `alpha_{t-1}`, `alpha_t` y `alpha_{t+1}`.

---

## 7. Bloques diagonales de `Omega`

Definamos, para abreviar,

\[
M=Z^\top D^{-1}Z.
\]

### Primer bloque diagonal

Para el primer estado,

\[
\boxed{
\Omega_{11}=Z^\top D^{-1}Z+T^\top Q^{-1}T+Q_1^{-1}
}
\]

O usando `M`,

\[
\Omega_{11}=M+T^\top Q^{-1}T+Q_1^{-1}.
\]

Los tres términos corresponden a:

1. La observación `y_1`.
2. La transición de `alpha_1` hacia `alpha_2`.
3. La distribución inicial de `alpha_1`.

### Bloques diagonales interiores

Para `t=2,...,n-1`,

\[
\boxed{
\Omega_{tt}=Z^\top D^{-1}Z+Q^{-1}+T^\top Q^{-1}T
}
\]

O equivalentemente,

\[
\Omega_{tt}=M+Q^{-1}+T^\top Q^{-1}T.
\]

Aquí aparecen:

1. La observación en el período `t`.
2. La transición que llega desde `alpha_{t-1}`.
3. La transición que sale hacia `alpha_{t+1}`.

### Último bloque diagonal

Para el último estado,

\[
\boxed{
\Omega_{nn}=Z^\top D^{-1}Z+Q^{-1}
}
\]

O equivalentemente,

\[
\Omega_{nn}=M+Q^{-1}.
\]

No aparece `T^T Q^{-1} T` porque no hay una transición posterior hacia `alpha_{n+1}`.

---

## 8. Bloques fuera de la diagonal

Los bloques que conectan períodos consecutivos son

\[
\boxed{
\Omega_{t,t+1}=-T^\top Q^{-1}
}
\]

\[
\boxed{
\Omega_{t+1,t}=-Q^{-1}T
}
\]

para `t=1,...,n-1`.

Como `Omega` es simétrica,

\[
\Omega_{t+1,t}=\Omega_{t,t+1}^\top.
\]

Todos los demás bloques son cero:

\[
\boxed{
\Omega_{ts}=0\quad\text{si }|t-s|>1.
}
\]

---

## 9. El vector `c`

El vector `c` se apila como

\[
c=
\begin{pmatrix}
 c_1\\
 c_2\\
\vdots\\
 c_n
\end{pmatrix},
\qquad
c_t\in\mathbb{R}^{m}.
\]

La contribución de las observaciones en cada período es

\[
Z^\top D^{-1}y_t.
\]

Por ello, para el primer período,

\[
\boxed{
 c_1=Z^\top D^{-1}y_1+Q_1^{-1}a
}
\]

El segundo término proviene de la distribución inicial

\[
\alpha_1\sim N(a,Q_1).
\]

Para `t=2,...,n`,

\[
\boxed{
 c_t=Z^\top D^{-1}y_t.
}
\]

Si `a=0`, entonces la fórmula se simplifica a

\[
\boxed{
 c_t=Z^\top D^{-1}y_t,
\qquad t=1,\ldots,n.
}
\]

---

## 10. Caso con intercepto en la transición

Si el modelo fuera

\[
\alpha_{t+1}=d+T\alpha_t+v_t,
\]

entonces el vector `c` recibiría términos adicionales:

\[
 c_1=Z^\top D^{-1}y_1+Q_1^{-1}a-T^\top Q^{-1}d,
\]

\[
 c_t=Z^\top D^{-1}y_t+Q^{-1}d-T^\top Q^{-1}d,
\qquad t=2,\ldots,n-1,
\]

\[
 c_n=Z^\top D^{-1}y_n+Q^{-1}d.
\]

En el modelo factorial dinámico básico de McCausland et al. se toma `d=0`, por lo que estos términos desaparecen.

---

## 11. Recuperación de las innovaciones después del suavizamiento

Una vez obtenida una trayectoria suavizada o simulada

\[
\alpha_1,\alpha_2,\ldots,\alpha_n,
\]

el shock del estado se recupera mediante

\[
\boxed{
 v_t=\alpha_{t+1}-T\alpha_t,
\qquad t=1,\ldots,n-1.
}
\]

El error de observación se recupera mediante

\[
\boxed{
 u_t=y_t-Z\alpha_t,
\qquad t=1,\ldots,n.
}
\]

Estos son los valores realizados de los errores una vez estimados los estados.

---

## 12. Ejemplo de dimensiones

Supongamos que hay:

- `p=10` series observadas;
- `m=3` factores;
- `n=100` períodos.

Entonces,

\[
y_t:10\times1,
\qquad
\alpha_t:3\times1,
\]

\[
Z:10\times3,
\qquad
D:10\times10,
\]

\[
T,Q,Q_1:3\times3.
\]

Cada bloque de `Omega` es `3 x 3`, pero la matriz completa es

\[
\Omega:300\times300.
\]

El vector completo es

\[
c:300\times1.
\]

La matriz `Omega` tiene la forma

\[
\Omega=
\begin{pmatrix}
3\times3&3\times3&0&\cdots&0\\
3\times3&3\times3&3\times3&\cdots&0\\
0&3\times3&3\times3&\cdots&0\\
\vdots&\vdots&\vdots&\ddots&3\times3\\
0&0&0&3\times3&3\times3
\end{pmatrix}.
\]

---

## 13. Relación con el artículo

En la aplicación de McCausland et al., esta estructura se usa para el suavizamiento eficiente de los estados latentes. Los autores comparan métodos basados en el filtro de Kalman con métodos basados en la matriz de precisión, incluyendo el algoritmo CFA y el algoritmo MMP.[web:97]

La ventaja de escribir el modelo mediante `Omega` y `c` es que no es necesario construir ni invertir directamente una matriz de covarianza densa de dimensión `nm x nm`. La estructura tridiagonal por bloques permite utilizar factorizaciones y sustituciones eficientes.

---

## 14. Fórmulas principales reunidas

Para el modelo

\[
y_t=Z\alpha_t+u_t,
\qquad u_t\sim N(0,D),
\]

\[
\alpha_1\sim N(a,Q_1),
\]

\[
\alpha_{t+1}=T\alpha_t+v_t,
\qquad v_t\sim N(0,Q),
\]

las expresiones principales son

\[
A_t=
\begin{pmatrix}
D^{-1}&0\\
0&Q^{-1}
\end{pmatrix},
\]

\[
\Omega_{11}=Z^\top D^{-1}Z+T^\top Q^{-1}T+Q_1^{-1},
\]

\[
\Omega_{tt}=Z^\top D^{-1}Z+Q^{-1}+T^\top Q^{-1}T,
\qquad t=2,\ldots,n-1,
\]

\[
\Omega_{nn}=Z^\top D^{-1}Z+Q^{-1},
\]

\[
\Omega_{t,t+1}=-T^\top Q^{-1},
\]

\[
\Omega_{t+1,t}=-Q^{-1}T,
\]

\[
c_1=Z^\top D^{-1}y_1+Q_1^{-1}a,
\]

\[
c_t=Z^\top D^{-1}y_t,
\qquad t=2,\ldots,n,
\]

\[
\alpha\mid y\sim N(\Omega^{-1}c,\Omega^{-1}),
\]

\[
E(\alpha\mid y)=\Omega^{-1}c.
\]
