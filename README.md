https://github.com/Sofia-Cadavid-Rueda/st-2026-2-tarea1-herramientas-Cadavid-Rueda-Sofia

# Tarea 1: Caja de herramientas de pronóstico

## Qué contiene el repositorio

| Archivo | Contenido | Dependencias |
|---|---|---|
| `R/00-lectura.R` | `leer_serie()` | R base |
| `R/01-graficos.R` | `graficar_serie()`, `correlograma()` | ggplot2, stats |
| `R/02-metodos.R` | los ocho métodos de pronóstico y `optimizar()` | tibble, stats |
| `R/03-evaluacion.R` | `medidas()`, `ljung_box()`, `jarque_bera()`, `durbin_watson()`, `validar_errores()` | stats |
| `ejemplos/ejemplos.R` | carga `R/` y corre los ocho ejemplos en orden | todo lo anterior, ggplot2, patchwork |
| `informe/informe.qmd` | análisis completo de los ocho ejemplos | dplyr, tidyr, purrr, tibble, ggplot2, patchwork |

## Cómo se corre

```r
source("ejemplos/ejemplos.R")
```

Tiempo de ejecución observado: **0.03 s** (medido con `system.time()`).

## Cómo se usan las funciones

```r
source("R/00-lectura.R")
source("R/02-metodos.R")

datos  <- leer_serie(AirPassengers, fuente = "AirPassengers, paquete datasets de R", unidad = "pasajeros (miles)")
ajuste <- ajustar_media(datos$y)
ajuste$pronosticar(6)
```

## Convenciones que fijan los números

- ACF calculada con divisor único `T` (Definición 2.5 de las notas de clase), no `T-h`.
- Media simple: calentamiento `Ŷ₂ = Y₁`; recursión `Ŷ_{t+1} = t⁻¹Σ Yᵢ` (no la media de toda la muestra).
- Media móvil: calentamiento de `k` períodos.
- SES: calentamiento `Ŷ₂ = Y₁`; forma de corrección de error `Ŷ_{t+1} = Ŷ_t + αe_t`, verificada contra la forma de pesos geométricos.
- Doble media móvil: calentamiento de `2k-1` períodos.
- Holt lineal: `L₁ = Y₁`, `T̂₁ = 0`.
- Tendencias: estimadas por ecuaciones normales; error estándar robusto HAC con núcleo de Bartlett y `⌊4(T/100)^(2/9)⌋` rezagos.
- Banda del correlograma: `±1.96/√n`, con `n` = observaciones de la serie o número de errores, según si se aplica a la serie o a los residuos de un método.

## Resumen de resultados

| Ejemplo | Serie | Parámetros | MASE método | MASE referente |
|---|---|---|---|---|
| Media simple | discoveries | — | 0.919 | 1.136 |
| Media móvil | discoveries | k = 9 | 0.588 | 1.136 |
| SES | LakeHuron | α = 0.98 | 2.141 | 2.164 |
| Doble media móvil | Nile | k = 12 | 0.951 | 0.835 |
| Tendencia lineal | austres | — | 3.974 | 6.091 |
| Tendencia cuadrática | co2 | — | 1.904 | 1.916 |
| Tendencia exponencial | JohnsonJohnson | α=—, β=— | 4.225 | 9.458 |
| Holt lineal | austres | α = 0.95, β = 0.80 | 1.484 | 6.091 |

## Declaración de uso de IA

Se usó Claude (Anthropic) como asistente durante el desarrollo de esta tarea.

Se pidió: explicación conceptual de cada método antes de implementarlo, plantillas de código para los ocho métodos y las funciones de evaluación, ayuda para depurar errores de ejecución y de renderizado en Quarto (paneles con pestañas anidadas, autoimpresión de expresiones sueltas en los chunks), y la estructura general del informe y del README.

Se recibió: explicaciones de cada fórmula contrastadas contra las notas de clase, código con verificación numérica incorporada contra funciones de R (`lm()`, `acf()`, `Box.test()`, pesos geométricos), y plantillas de chunks de código para el informe y para `ejemplos.R`.

Se verificó por cuenta propia: cada función se corrió sobre datos reales y sus salidas se compararon contra los ejemplos numéricos de las notas de clase antes de aceptarlas (incluyendo un error real detectado por la estudiante en la implementación de `medidas()`, corregido durante el desarrollo); los análisis en prosa de cada ejemplo (lectura del patrón, interpretación de las pruebas de hipótesis, conclusiones y recomendaciones) fueron escritos por la estudiante.