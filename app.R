# Simulaciones de procesos estadísticos aleatorios
# - Ley de los grandes números
# - Aproximación de la distribución binomial mediante la Poisson
# - Teorema central del límite
# - Intervalos de confianza para la media y para proporciones
# - Contrastes de hipótesis para la comparación de dos medias o dos proporciones

library(shiny)
library(bslib)
library(ggplot2)

# Carpeta de imágenes (logos) accesible desde el navegador en la ruta "img/".
addResourcePath("img", "img")

# Paleta y tema común de los gráficos -----------------------------------------

# Colores de los gráficos para el tema claro y el oscuro. Las series son los
# mismos tonos ajustados a cada fondo (azul, naranja, aguamarina, ...).
modos <- list(
  claro = list(
    series = c(
      "#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4",
      "#008300", "#6250d6", "#e34948"
    ),
    texto = "#212529", secundario = "#52514e", rejilla = "#e6e6e3",
    superficie = "#ffffff"
  ),
  oscuro = list(
    series = c(
      "#3987e5", "#d95926", "#199e70", "#c98500", "#d55181",
      "#008300", "#9085e9", "#e66767"
    ),
    texto = "#dee2e6", secundario = "#c3c2b7", rejilla = "#3a3c3f",
    superficie = "#212529"
  )
)

# El fondo del gráfico es el mismo que el de las tarjetas en cada modo.
tema <- function(m) {
  theme_minimal(base_size = 14) +
    theme(
      plot.background = element_rect(fill = m$superficie, colour = NA),
      text = element_text(colour = m$texto),
      axis.text = element_text(colour = m$secundario),
      axis.title = element_text(colour = m$secundario),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = m$rejilla, linewidth = 0.3),
      legend.position = "top",
      plot.title = element_text(face = "bold")
    )
}

fmt <- function(x, d = 4) formatC(x, digits = d, format = "f", decimal.mark = ",")
pct <- function(x, d = 1) paste0(fmt(100 * x, d), " %")

# Distribuciones de población para el teorema central del límite --------------

poblacion <- function(input) {
  switch(input$tcl_dist,
    uniforme = {
      a <- input$unif_a
      b <- max(input$unif_b, a + 0.01)
      list(
        nombre = sprintf("U(%s, %s)", a, b), discreta = FALSE,
        r = function(k) runif(k, a, b),
        d = function(x) dunif(x, a, b),
        rango = c(a - 0.1 * (b - a), b + 0.1 * (b - a)),
        media = (a + b) / 2, sd = (b - a) / sqrt(12)
      )
    },
    binomial = {
      n <- input$binom_n
      p <- input$binom_p
      list(
        nombre = sprintf("B(%d, %s)", n, p), discreta = TRUE,
        r = function(k) rbinom(k, n, p),
        d = function(x) dbinom(x, n, p),
        rango = c(0, n),
        media = n * p, sd = sqrt(n * p * (1 - p))
      )
    },
    poisson = {
      l <- input$pois_lambda
      list(
        nombre = sprintf("P(%s)", l), discreta = TRUE,
        r = function(k) rpois(k, l),
        d = function(x) dpois(x, l),
        rango = c(0, qpois(0.9995, l)),
        media = l, sd = sqrt(l)
      )
    },
    normal = {
      m <- input$norm_mu
      s <- input$norm_sigma
      list(
        nombre = sprintf("N(%s, %s)", m, s), discreta = FALSE,
        r = function(k) rnorm(k, m, s),
        d = function(x) dnorm(x, m, s),
        rango = c(m - 4 * s, m + 4 * s),
        media = m, sd = s
      )
    },
    exponencial = {
      l <- input$exp_lambda
      list(
        nombre = sprintf("Exp(%s)", l), discreta = FALSE,
        r = function(k) rexp(k, l),
        d = function(x) dexp(x, l),
        rango = c(0, qexp(0.999, l)),
        media = 1 / l, sd = 1 / l
      )
    },
    chicuadrado = {
      g <- input$chi_gl
      list(
        nombre = sprintf("χ²(%d)", g), discreta = FALSE,
        r = function(k) rchisq(k, g),
        d = function(x) dchisq(x, g),
        rango = c(0, qchisq(0.999, g)),
        media = g, sd = sqrt(2 * g)
      )
    }
  )
}

# Contrastes de comparación de dos poblaciones ---------------------------------

# Los parámetros de un contraste son una lista con el tipo ("media" o
# "proporcion"), los tamaños muestrales n1 y n2, el nivel de significación alfa,
# la hipótesis alternativa alt ("bilateral", "mayor" o "menor") y los parámetros
# de cada población: mu1, mu2, sd1, sd2 y var ("welch" o "iguales") en las
# medias, o p1 y p2 en las proporciones.

# Bajo la hipótesis nula la segunda población tiene la misma media (o
# proporción) que la primera; el resto de parámetros no cambia.
bajo_h0 <- function(s) {
  if (s$tipo == "media") s$mu2 <- s$mu1 else s$p2 <- s$p1
  s
}

diferencia <- function(s) if (s$tipo == "media") s$mu1 - s$mu2 else s$p1 - s$p2

# Grados de libertad de Welch-Satterthwaite a partir de las varianzas v1 y v2.
gl_welch <- function(v1, v2, n1, n2) {
  a <- v1 / n1
  b <- v2 / n2
  (a + b)^2 / (a^2 / (n1 - 1) + b^2 / (n2 - 1))
}

# Rechazo de H0 según la hipótesis alternativa, con q la función cuantil de la
# distribución del estadístico bajo H0.
rechazo <- function(t, alt, alfa, q) {
  switch(alt,
    bilateral = abs(t) > q(1 - alfa / 2),
    mayor = t > q(1 - alfa),
    menor = t < q(alfa)
  )
}

# Simula k contrastes y devuelve el estadístico y la decisión de cada uno. En
# las medias de poblaciones normales se simulan directamente la media y la
# cuasivarianza de cada muestra (x̄ ~ N(μ, σ/√n) y (n − 1)s²/σ² ~ χ²(n − 1),
# independientes), lo que equivale a generar las muestras y es mucho más rápido.
contraste <- function(s, k) {
  n1 <- s$n1
  n2 <- s$n2
  if (s$tipo == "media") {
    m1 <- rnorm(k, s$mu1, s$sd1 / sqrt(n1))
    m2 <- rnorm(k, s$mu2, s$sd2 / sqrt(n2))
    v1 <- s$sd1^2 * rchisq(k, n1 - 1) / (n1 - 1)
    v2 <- s$sd2^2 * rchisq(k, n2 - 1) / (n2 - 1)
    if (s$var == "iguales") {
      gl <- n1 + n2 - 2
      se <- sqrt(((n1 - 1) * v1 + (n2 - 1) * v2) / gl * (1 / n1 + 1 / n2))
    } else {
      gl <- gl_welch(v1, v2, n1, n2)
      se <- sqrt(v1 / n1 + v2 / n2)
    }
    t <- (m1 - m2) / se
    q <- function(p) qt(p, gl)
  } else {
    x1 <- rbinom(k, n1, s$p1)
    x2 <- rbinom(k, n2, s$p2)
    p <- (x1 + x2) / (n1 + n2)
    se <- sqrt(p * (1 - p) * (1 / n1 + 1 / n2))
    # Si todas las observaciones son éxitos (o fracasos) no hay diferencia.
    t <- ifelse(se > 0, (x1 / n1 - x2 / n2) / se, 0)
    q <- qnorm
  }
  list(t = t, rechaza = rechazo(t, s$alt, s$alfa, q))
}

# Valores críticos del estadístico. En el contraste de Welch los grados de
# libertad dependen de cada muestra; se usan los de las varianzas poblacionales.
criticos <- function(s) {
  q <- if (s$tipo == "proporcion") {
    qnorm
  } else if (s$var == "iguales") {
    function(p) qt(p, s$n1 + s$n2 - 2)
  } else {
    function(p) qt(p, gl_welch(s$sd1^2, s$sd2^2, s$n1, s$n2))
  }
  switch(s$alt,
    bilateral = c(-1, 1) * q(1 - s$alfa / 2),
    mayor = q(1 - s$alfa),
    menor = q(s$alfa)
  )
}

# Distribución teórica del estadístico de contraste: lista con la densidad d,
# la función de distribución p y la función cuantil q.
# - Medias: el estadístico es (aproximadamente) escala · t no central con gl
#   grados de libertad y parámetro de no centralidad (μ₁ − μ₂)/√(σ₁²/n₁ + σ₂²/n₂).
#   Con la t de Student y σ₁ = σ₂ es exacta (escala = 1, gl = n₁ + n₂ − 2). Con
#   σ₁ ≠ σ₂ la varianza ponderada no estima el error típico real, lo que da la
#   escala, y sus grados de libertad se aproximan por Satterthwaite. En la t de
#   Welch se usan los grados de libertad de Welch-Satterthwaite.
# - Proporciones: aproximación normal de la diferencia de proporciones
#   muestrales, con el error típico ponderado se0 en el denominador del
#   estadístico y el error típico real se1.
distribucion <- function(s) {
  if (s$tipo == "media") {
    se <- sqrt(s$sd1^2 / s$n1 + s$sd2^2 / s$n2)
    ncp <- diferencia(s) / se
    if (s$var == "iguales") {
      a <- (s$n1 - 1) * s$sd1^2
      b <- (s$n2 - 1) * s$sd2^2
      gl <- (a + b)^2 / ((s$n1 - 1) * s$sd1^4 + (s$n2 - 1) * s$sd2^4)
      escala <- se / sqrt((a + b) / (s$n1 + s$n2 - 2) * (1 / s$n1 + 1 / s$n2))
    } else {
      gl <- gl_welch(s$sd1^2, s$sd2^2, s$n1, s$n2)
      escala <- 1
    }
    # La t no central avisa de pérdidas de precisión irrelevantes aquí.
    list(
      d = function(x) suppressWarnings(dt(x / escala, gl, ncp)) / escala,
      p = function(x) suppressWarnings(pt(x / escala, gl, ncp)),
      q = function(p) suppressWarnings(qt(p, gl, ncp)) * escala
    )
  } else {
    p <- (s$n1 * s$p1 + s$n2 * s$p2) / (s$n1 + s$n2)
    se0 <- sqrt(p * (1 - p) * (1 / s$n1 + 1 / s$n2))
    se1 <- sqrt(s$p1 * (1 - s$p1) / s$n1 + s$p2 * (1 - s$p2) / s$n2)
    media <- diferencia(s) / se0
    sd <- se1 / se0
    list(
      d = function(x) dnorm(x, media, sd),
      p = function(x) pnorm(x, media, sd),
      q = function(p) qnorm(p, media, sd)
    )
  }
}

# Probabilidad teórica de rechazar H0 con los parámetros dados: error de tipo I
# si H0 es cierta y potencia si no lo es.
prob_rechazo <- function(s) {
  f <- distribucion(s)$p
  crit <- criticos(s)
  switch(s$alt,
    bilateral = f(crit[1]) + 1 - f(crit[2]),
    mayor = 1 - f(crit),
    menor = f(crit)
  )
}

# Tamaño del efecto: d de Cohen con la desviación típica combinada
# √((σ₁² + σ₂²)/2) para las medias y h de Cohen para las proporciones.
tamano_efecto <- function(s) {
  if (s$tipo == "media") {
    diferencia(s) / sqrt((s$sd1^2 + s$sd2^2) / 2)
  } else {
    2 * asin(sqrt(s$p1)) - 2 * asin(sqrt(s$p2))
  }
}

magnitud_efecto <- function(e) {
  e <- abs(e)
  if (e < 0.2) "Insignificante" else if (e < 0.5) "Pequeño" else if (e < 0.8) "Mediano" else "Grande"
}

# Interfaz ---------------------------------------------------------------------

boton_simular <- function(id) {
  actionButton(id, "Simular de nuevo", icon = icon("dice"), class = "btn-primary")
}

# Fila de recuadros de valores: en pantallas pequeñas van de dos en dos (el
# último ocupa toda la fila si quedan impares) y en grandes, todos en una fila.
fila_valores <- function(...) {
  k <- length(list(...))
  movil <- rep(6, k)
  if (k %% 2 == 1) movil[k] <- 12
  layout_columns(
    fill = FALSE,
    class = "fila-valores",
    col_widths = breakpoints(xs = movil, lg = rep(12 / k, k)),
    ...
  )
}

# Recuadros de valores más compactos que los de bslib, y aún más en pantallas
# pequeñas (< 576 px). En móvil bslib fuerza todas las rejillas a una columna;
# en las filas de valores se restablecen las 12 columnas para que los
# recuadros vayan de dos en dos. Los selectores llevan la clase .fila-valores
# para ganar en especificidad al CSS de bslib, que se carga después.
estilo_movil <- "
.fila-valores .bslib-value-box > .card-body { padding: 0; }
.fila-valores .bslib-value-box .value-box-area { padding: 0.75rem 1rem; }
.fila-valores .bslib-value-box .value-box-title {
  font-size: 0.9rem; margin-bottom: 0.25rem;
}
.fila-valores .bslib-value-box .value-box-value {
  font-size: 1.5rem; margin-bottom: 0.15rem;
}
.fila-valores .bslib-value-box .value-box-area p { font-size: 0.8rem; }
@media (max-width: 575.98px) {
  .bslib-grid.fila-valores {
    grid-template-columns: repeat(12, minmax(0, 1fr)) !important;
    gap: 0.5rem;
  }
  .fila-valores .bslib-value-box { margin-bottom: 0; }
  .fila-valores .bslib-value-box .value-box-area { padding: 0.5rem 0.75rem; }
  .fila-valores .bslib-value-box .value-box-title {
    font-size: 0.75rem; line-height: 1.2; margin-bottom: 0.25rem;
  }
  .fila-valores .bslib-value-box .value-box-value { font-size: 1.1rem; }
}
"

ui_lgn <- nav_panel(
  "Ley de los grandes números",
  layout_sidebar(
    sidebar = sidebar(
      width = 300,
      title = "Parámetros",
      radioButtons("lgn_exp", "Experimento",
        choices = c(
          "Lanzar una moneda (proporción de caras)" = "moneda",
          "Lanzar un dado (media de los resultados)" = "dado"
        )
      ),
      conditionalPanel(
        "input.lgn_exp == 'moneda'",
        sliderInput("lgn_p", "Probabilidad de cara (p)", min = 0.01, max = 0.99, value = 0.5, step = 0.01)
      ),
      sliderInput("lgn_n", "Número de repeticiones", min = 10, max = 10000, value = 1000, step = 10),
      sliderInput("lgn_tray", "Número de secuencias simuladas", min = 1, max = 8, value = 3, step = 1),
      checkboxInput("lgn_log", "Escala logarítmica en el eje X", value = TRUE),
      boton_simular("lgn_sim")
    ),
    fila_valores(
      value_box("Media teórica (μ)", textOutput("lgn_mu"), theme = "primary"),
      value_box("Media simulada final (promedio de secuencias)", textOutput("lgn_media")),
      value_box("Error absoluto medio final", textOutput("lgn_error"))
    ),
    card(
      full_screen = TRUE,
      card_header("Evolución de la media muestral al aumentar el número de repeticiones"),
      plotOutput("lgn_plot", height = "450px")
    )
  )
)

ui_bp <- nav_panel(
  "Ley de los casos raros",
  layout_sidebar(
    sidebar = sidebar(
      width = 300,
      title = "Parámetros",
      sliderInput("bp_n", "Número de ensayos (n)", min = 1, max = 500, value = 50, step = 1),
      sliderInput("bp_p", "Probabilidad de éxito (p)", min = 0.001, max = 0.5, value = 0.05, step = 0.001),
      checkboxInput("bp_sim", "Mostrar también una simulación de la binomial", value = FALSE),
      conditionalPanel(
        "input.bp_sim",
        sliderInput("bp_m", "Número de valores simulados", min = 100, max = 10000, value = 1000, step = 100),
        boton_simular("bp_resim")
      ),
      helpText(
        "Regla práctica: la aproximación B(n, p) ≈ P(λ = np) es buena cuando",
        "n ≥ 30 y p ≤ 0,1 (n grande y p pequeño)."
      )
    ),
    fila_valores(
      value_box("λ = n·p", textOutput("bp_lambda"), theme = "primary"),
      value_box("Varianza binomial n·p·(1−p)", textOutput("bp_var")),
      value_box("Máxima diferencia de probabilidades", textOutput("bp_maxdif")),
      value_box("¿Se cumple la regla práctica?", textOutput("bp_regla"))
    ),
    card(
      full_screen = TRUE,
      card_header("Función de probabilidad de la binomial y de la Poisson"),
      plotOutput("bp_plot", height = "450px")
    )
  )
)

ui_tcl <- nav_panel(
  "Teorema central del límite",
  layout_sidebar(
    sidebar = sidebar(
      width = 300,
      title = "Parámetros",
      selectInput("tcl_dist", "Distribución de la población",
        choices = c(
          "Uniforme" = "uniforme", "Binomial" = "binomial", "Poisson" = "poisson",
          "Normal" = "normal", "Exponencial" = "exponencial", "Chi-cuadrado" = "chicuadrado"
        )
      ),
      conditionalPanel(
        "input.tcl_dist == 'uniforme'",
        numericInput("unif_a", "Mínimo (a)", value = 0),
        numericInput("unif_b", "Máximo (b)", value = 1)
      ),
      conditionalPanel(
        "input.tcl_dist == 'binomial'",
        sliderInput("binom_n", "Número de ensayos", min = 1, max = 50, value = 10, step = 1),
        sliderInput("binom_p", "Probabilidad de éxito", min = 0.01, max = 0.99, value = 0.2, step = 0.01)
      ),
      conditionalPanel(
        "input.tcl_dist == 'poisson'",
        sliderInput("pois_lambda", "λ", min = 0.1, max = 20, value = 2, step = 0.1)
      ),
      conditionalPanel(
        "input.tcl_dist == 'normal'",
        numericInput("norm_mu", "Media (μ)", value = 0),
        numericInput("norm_sigma", "Desviación típica (σ)", value = 1, min = 0.01)
      ),
      conditionalPanel(
        "input.tcl_dist == 'exponencial'",
        sliderInput("exp_lambda", "λ", min = 0.1, max = 5, value = 1, step = 0.1)
      ),
      conditionalPanel(
        "input.tcl_dist == 'chicuadrado'",
        sliderInput("chi_gl", "Grados de libertad", min = 1, max = 30, value = 2, step = 1)
      ),
      sliderInput("tcl_n", "Tamaño muestral (n)", min = 1, max = 200, value = 30, step = 1),
      sliderInput("tcl_m", "Número de muestras", min = 100, max = 10000, value = 2000, step = 100),
      checkboxInput("tcl_normal", "Superponer la normal N(μ, σ/√n)", value = TRUE),
      boton_simular("tcl_sim")
    ),
    fila_valores(
      value_box("Media teórica de X̄ (μ)", textOutput("tcl_mu"), theme = "primary"),
      value_box("Media de las medias simuladas", textOutput("tcl_media")),
      value_box("Desv. típica teórica de X̄ (σ/√n)", textOutput("tcl_sd"), theme = "primary"),
      value_box("Desv. típica de las medias simuladas", textOutput("tcl_sdsim"))
    ),
    layout_columns(
      col_widths = c(8, 4),
      card(
        full_screen = TRUE,
        card_header("Distribución de la media muestral X̄"),
        plotOutput("tcl_plot", height = "450px")
      ),
      card(
        full_screen = TRUE,
        card_header(textOutput("tcl_pob_titulo", inline = TRUE)),
        plotOutput("tcl_pob_plot", height = "450px")
      )
    )
  )
)

ui_ic <- nav_panel(
  "Intervalos de confianza",
  layout_sidebar(
    sidebar = sidebar(
      width = 300,
      title = "Parámetros",
      radioButtons("ic_tipo", "Parámetro a estimar",
        choices = c("Media de una población normal" = "media", "Proporción" = "proporcion")
      ),
      conditionalPanel(
        "input.ic_tipo == 'media'",
        numericInput("ic_mu", "Media de la población (μ)", value = 100),
        numericInput("ic_sigma", "Desviación típica de la población (σ)", value = 15, min = 0.01),
        radioButtons("ic_sigma_conocida", "Desviación típica poblacional",
          choices = c("Conocida (intervalo z)" = "z", "Desconocida (intervalo t)" = "t")
        )
      ),
      conditionalPanel(
        "input.ic_tipo == 'proporcion'",
        sliderInput("ic_p", "Proporción poblacional (p)", min = 0.01, max = 0.99, value = 0.3, step = 0.01),
        radioButtons("ic_metodo", "Método",
          choices = c("Wald (aproximación normal)" = "wald", "Wilson" = "wilson")
        )
      ),
      sliderInput("ic_n", "Tamaño muestral (n)", min = 2, max = 500, value = 30, step = 1),
      sliderInput("ic_conf", "Nivel de confianza", min = 0.80, max = 0.99, value = 0.95, step = 0.01),
      sliderInput("ic_k", "Número de intervalos", min = 10, max = 500, value = 100, step = 10),
      boton_simular("ic_sim")
    ),
    fila_valores(
      value_box("Nivel de confianza", textOutput("ic_nivel"), theme = "primary"),
      value_box("Intervalos que contienen el parámetro", textOutput("ic_cobertura")),
      value_box("Amplitud media de los intervalos", textOutput("ic_amplitud"))
    ),
    layout_columns(
      col_widths = c(7, 5),
      card(
        full_screen = TRUE,
        card_header("Intervalos de confianza obtenidos en cada muestra"),
        plotOutput("ic_plot", height = "500px")
      ),
      card(
        full_screen = TRUE,
        card_header("Proporción acumulada de intervalos que contienen el parámetro"),
        plotOutput("ic_cob_plot", height = "500px")
      )
    )
  )
)

ui_ct <- nav_panel(
  "Contrastes de hipótesis",
  layout_sidebar(
    sidebar = sidebar(
      width = 300,
      title = "Parámetros",
      radioButtons("ct_tipo", "Comparación",
        choices = c(
          "Medias de dos poblaciones normales" = "media",
          "Proporciones de dos poblaciones" = "proporcion"
        )
      ),
      conditionalPanel(
        "input.ct_tipo == 'media'",
        numericInput("ct_mu1", "Media de la población 1 (μ₁)", value = 100),
        numericInput("ct_mu2", "Media de la población 2 (μ₂)", value = 105),
        sliderInput("ct_sd1", "Desviación típica de la población 1 (σ₁)", min = 1, max = 50, value = 15, step = 0.5),
        sliderInput("ct_sd2", "Desviación típica de la población 2 (σ₂)", min = 1, max = 50, value = 15, step = 0.5),
        radioButtons("ct_var", "Varianzas poblacionales (desconocidas)",
          choices = c(
            "Distintas (t de Welch)" = "welch",
            "Iguales (t de Student con varianza ponderada)" = "iguales"
          )
        )
      ),
      conditionalPanel(
        "input.ct_tipo == 'proporcion'",
        sliderInput("ct_p1", "Proporción de la población 1 (p₁)", min = 0.01, max = 0.99, value = 0.3, step = 0.01),
        sliderInput("ct_p2", "Proporción de la población 2 (p₂)", min = 0.01, max = 0.99, value = 0.4, step = 0.01)
      ),
      sliderInput("ct_n1", "Tamaño de la muestra 1 (n₁)", min = 2, max = 500, value = 30, step = 1),
      sliderInput("ct_n2", "Tamaño de la muestra 2 (n₂)", min = 2, max = 500, value = 30, step = 1),
      selectInput("ct_alt", "Hipótesis alternativa (H₁)",
        choices = c(
          "Bilateral: población 1 ≠ población 2" = "bilateral",
          "Unilateral: población 1 > población 2" = "mayor",
          "Unilateral: población 1 < población 2" = "menor"
        )
      ),
      sliderInput("ct_alfa", "Nivel de significación (α)", min = 0.01, max = 0.20, value = 0.05, step = 0.01),
      sliderInput("ct_k", "Número de contrastes simulados", min = 100, max = 10000, value = 2000, step = 100),
      radioButtons("ct_distr", "Distribución del estadístico",
        choices = c("Simulada (histograma)" = "sim", "Teórica (curva)" = "teo", "Ambas" = "ambas")
      ),
      radioButtons("ct_eje", "En el gráfico de errores, variar",
        choices = c("El tamaño muestral" = "n", "La variabilidad de las poblaciones" = "var")
      ),
      checkboxInput("ct_teorica", "Superponer los valores teóricos", value = TRUE),
      boton_simular("ct_sim"),
      helpText(
        "Para estimar el error de tipo I se simulan contrastes con H₀ cierta",
        "(μ₂ = μ₁ o p₂ = p₁, manteniendo el resto de parámetros); para el error",
        "de tipo II y la potencia, con los parámetros elegidos."
      ),
      helpText(
        "Las distribuciones teóricas son exactas para la t de Student con σ₁ = σ₂.",
        "En los demás casos son aproximaciones: t no central con grados de",
        "libertad de Satterthwaite para las medias y normal para las proporciones."
      )
    ),
    fila_valores(
      value_box(
        textOutput("ct_efecto_titulo"), textOutput("ct_efecto"),
        p(textOutput("ct_efecto_magnitud", inline = TRUE)),
        theme = "primary"
      ),
      value_box("Error de tipo I simulado", textOutput("ct_tipo1"), p(textOutput("ct_tipo1_teo", inline = TRUE))),
      value_box("Error de tipo II simulado (β)", textOutput("ct_tipo2"), p(textOutput("ct_tipo2_teo", inline = TRUE))),
      value_box("Potencia simulada (1 − β)", textOutput("ct_potencia"), p(textOutput("ct_potencia_teo", inline = TRUE)))
    ),
    layout_columns(
      col_widths = c(6, 6),
      card(
        full_screen = TRUE,
        card_header("Distribución del estadístico de contraste bajo H₀ y bajo H₁"),
        plotOutput("ct_plot", height = "450px")
      ),
      card(
        full_screen = TRUE,
        card_header(textOutput("ct_curva_titulo", inline = TRUE)),
        plotOutput("ct_curva_plot", height = "450px")
      )
    )
  )
)

ui <- page_navbar(
  title = tags$span(
    class = "d-inline-flex align-items-center gap-2",
    tags$img(src = "img/logos/sticker.png", height = "36px", alt = ""),
    "Simulaciones estadísticas"
  ),
  window_title = "Simulaciones estadísticas",
  theme = bs_theme(version = 5, primary = modos$claro$series[1]) |> bs_add_rules(estilo_movil),
  ui_lgn,
  ui_bp,
  ui_tcl,
  ui_ic,
  ui_ct,
  nav_spacer(),
  nav_item(input_dark_mode(id = "modo", mode = "light"))
)

# Servidor ---------------------------------------------------------------------

server <- function(input, output, session) {
  m <- reactive(if (identical(input$modo, "dark")) modos$oscuro else modos$claro)

  # Ley de los grandes números -------------------------------------------------

  lgn <- reactive({
    input$lgn_sim
    n <- input$lgn_n
    if (input$lgn_exp == "moneda") {
      p <- input$lgn_p
      gen <- function() rbinom(n, 1, p)
      mu <- p
      etiqueta <- "Proporción de caras"
    } else {
      gen <- function() sample(1:6, n, replace = TRUE)
      mu <- 3.5
      etiqueta <- "Media de los resultados"
    }
    datos <- do.call(rbind, lapply(seq_len(input$lgn_tray), function(j) {
      x <- gen()
      data.frame(secuencia = factor(j), n = seq_len(n), media = cumsum(x) / seq_len(n))
    }))
    list(datos = datos, mu = mu, etiqueta = etiqueta)
  })

  finales <- reactive({
    d <- lgn()$datos
    d[d$n == max(d$n), "media"]
  })

  output$lgn_mu <- renderText(fmt(lgn()$mu))
  output$lgn_media <- renderText(fmt(mean(finales())))
  output$lgn_error <- renderText(fmt(mean(abs(finales() - lgn()$mu))))

  output$lgn_plot <- renderPlot({
    s <- lgn()
    m <- m()
    g <- ggplot(s$datos, aes(n, media, colour = secuencia)) +
      geom_hline(yintercept = s$mu, linetype = "dashed", linewidth = 0.6, colour = m$secundario) +
      geom_line(linewidth = 0.7) +
      scale_colour_manual(values = m$series, name = "Secuencia") +
      labs(
        x = "Número de repeticiones", y = s$etiqueta,
        subtitle = paste("Línea discontinua: valor teórico μ =", fmt(s$mu, 2))
      ) +
      tema(m)
    if (input$lgn_log) g <- g + scale_x_log10()
    if (input$lgn_tray == 1) g <- g + theme(legend.position = "none")
    g
  })

  # Aproximación binomial → Poisson --------------------------------------------

  bp <- reactive({
    n <- input$bp_n
    p <- input$bp_p
    l <- n * p
    xmax <- max(qbinom(0.9999, n, p), qpois(0.9999, l), 5)
    x <- 0:xmax
    list(
      n = n, p = p, l = l, x = x,
      binom = dbinom(x, n, p), pois = dpois(x, l)
    )
  })

  bp_muestra <- reactive({
    input$bp_resim
    rbinom(input$bp_m, input$bp_n, input$bp_p)
  })

  output$bp_lambda <- renderText(fmt(bp()$l))
  output$bp_var <- renderText(paste0(fmt(bp()$n * bp()$p * (1 - bp()$p)), "  (Poisson: ", fmt(bp()$l), ")"))
  output$bp_maxdif <- renderText(fmt(max(abs(bp()$binom - bp()$pois))))
  output$bp_regla <- renderText(if (bp()$n >= 30 && bp()$p <= 0.1) "Sí" else "No")

  output$bp_plot <- renderPlot({
    s <- bp()
    m <- m()
    d <- rbind(
      data.frame(x = s$x, prob = s$binom, dist = sprintf("Binomial B(%d, %s)", s$n, s$p)),
      data.frame(x = s$x, prob = s$pois, dist = sprintf("Poisson P(%s)", fmt(s$l, 3)))
    )
    d$dist <- factor(d$dist, levels = unique(d$dist))
    colores <- setNames(m$series[1:2], levels(d$dist))
    g <- ggplot(d, aes(x, prob, fill = dist)) +
      geom_col(position = position_dodge(width = 0.85), width = 0.8) +
      scale_fill_manual(values = colores, name = NULL) +
      labs(x = "Número de éxitos (x)", y = "P(X = x)") +
      tema(m)
    if (input$bp_sim) {
      muestra <- bp_muestra()
      f <- as.data.frame(table(factor(muestra, levels = s$x)), stringsAsFactors = FALSE)
      f$x <- as.numeric(f$Var1)
      f$prob <- f$Freq / length(muestra)
      f <- f[f$Freq > 0, ]
      g <- g +
        geom_point(
          data = f, aes(x, prob, colour = "Frecuencias relativas simuladas"),
          inherit.aes = FALSE, size = 3, shape = 21, fill = m$series[3], stroke = 0.8
        ) +
        scale_colour_manual(values = c("Frecuencias relativas simuladas" = m$superficie), name = NULL) +
        guides(colour = guide_legend(override.aes = list(fill = m$series[3], colour = m$secundario)))
    }
    if (length(s$x) <= 30) g <- g + scale_x_continuous(breaks = s$x)
    g
  })

  # Teorema central del límite -------------------------------------------------

  pob <- reactive(poblacion(input))

  tcl <- reactive({
    input$tcl_sim
    p <- pob()
    n <- input$tcl_n
    m <- input$tcl_m
    medias <- rowMeans(matrix(p$r(n * m), nrow = m, ncol = n))
    list(medias = medias, n = n, mu = p$media, sd = p$sd / sqrt(n), discreta = p$discreta)
  })

  output$tcl_mu <- renderText(fmt(tcl()$mu))
  output$tcl_media <- renderText(fmt(mean(tcl()$medias)))
  output$tcl_sd <- renderText(fmt(tcl()$sd))
  output$tcl_sdsim <- renderText(fmt(sd(tcl()$medias)))

  output$tcl_plot <- renderPlot({
    s <- tcl()
    m <- m()
    d <- data.frame(media = s$medias)
    # Ancho de clase de Freedman-Diaconis; en poblaciones discretas la media
    # sólo toma valores múltiplos de 1/n, así que se ajusta a ese retículo.
    ancho <- 2 * IQR(s$medias) / length(s$medias)^(1 / 3)
    if (!is.finite(ancho) || ancho <= 0) ancho <- max(s$sd / 10, 1e-6)
    if (s$discreta) {
      paso <- 1 / s$n
      ancho <- max(1, round(ancho / paso)) * paso
      limite <- paso / 2
    } else {
      limite <- 0
    }
    g <- ggplot(d, aes(media)) +
      geom_histogram(aes(y = after_stat(density)),
        binwidth = ancho, boundary = limite,
        fill = m$series[1], colour = m$superficie, linewidth = 0.3
      ) +
      geom_vline(xintercept = s$mu, linetype = "dashed", colour = m$secundario, linewidth = 0.6) +
      labs(x = "Media muestral (x̄)", y = "Densidad") +
      tema(m)
    if (input$tcl_normal) {
      g <- g + stat_function(
        fun = dnorm, args = list(mean = s$mu, sd = s$sd),
        colour = m$series[2], linewidth = 1, n = 500
      )
    }
    g
  })

  output$tcl_pob_titulo <- renderText(paste("Población:", pob()$nombre))

  output$tcl_pob_plot <- renderPlot({
    p <- pob()
    m <- m()
    if (p$discreta) {
      x <- seq(floor(p$rango[1]), ceiling(p$rango[2]))
      g <- ggplot(data.frame(x = x, y = p$d(x)), aes(x, y)) +
        geom_col(fill = m$series[3], width = 0.7) +
        labs(x = "x", y = "P(X = x)")
    } else {
      x <- seq(p$rango[1], p$rango[2], length.out = 500)
      g <- ggplot(data.frame(x = x, y = p$d(x)), aes(x, y)) +
        geom_area(fill = m$series[3], alpha = 0.3) +
        geom_line(colour = m$series[3], linewidth = 1) +
        labs(x = "x", y = "Densidad")
    }
    g +
      geom_vline(xintercept = p$media, linetype = "dashed", colour = m$secundario, linewidth = 0.6) +
      labs(subtitle = sprintf("μ = %s   σ = %s", fmt(p$media), fmt(p$sd))) +
      tema(m)
  })

  # Intervalos de confianza ----------------------------------------------------

  ic <- reactive({
    input$ic_sim
    k <- input$ic_k
    n <- input$ic_n
    alfa <- 1 - input$ic_conf
    z <- qnorm(1 - alfa / 2)
    if (input$ic_tipo == "media") {
      param <- input$ic_mu
      sigma <- max(input$ic_sigma, 0.01)
      x <- matrix(rnorm(n * k, param, sigma), nrow = k, ncol = n)
      est <- rowMeans(x)
      if (input$ic_sigma_conocida == "z") {
        error <- z * sigma / sqrt(n)
      } else {
        error <- qt(1 - alfa / 2, n - 1) * apply(x, 1, sd) / sqrt(n)
      }
      li <- est - error
      ls <- est + error
      simbolo <- "μ"
    } else {
      param <- input$ic_p
      exitos <- rbinom(k, n, param)
      est <- exitos / n
      if (input$ic_metodo == "wald") {
        error <- z * sqrt(est * (1 - est) / n)
        li <- est - error
        ls <- est + error
      } else {
        centro <- (exitos + z^2 / 2) / (n + z^2)
        error <- z * sqrt(n) / (n + z^2) * sqrt(est * (1 - est) + z^2 / (4 * n))
        li <- centro - error
        ls <- centro + error
      }
      li <- pmax(li, 0)
      ls <- pmin(ls, 1)
      simbolo <- "p"
    }
    datos <- data.frame(muestra = seq_len(k), est = est, li = li, ls = ls)
    datos$contiene <- factor(
      ifelse(li <= param & param <= ls, "Contiene el parámetro", "No contiene el parámetro"),
      levels = c("Contiene el parámetro", "No contiene el parámetro")
    )
    list(datos = datos, param = param, simbolo = simbolo, conf = input$ic_conf)
  })

  output$ic_nivel <- renderText(paste0(fmt(100 * ic()$conf, 0), " %"))
  output$ic_cobertura <- renderText({
    d <- ic()$datos
    aciertos <- sum(d$contiene == "Contiene el parámetro")
    sprintf("%d de %d (%s %%)", aciertos, nrow(d), fmt(100 * aciertos / nrow(d), 1))
  })
  output$ic_amplitud <- renderText(fmt(mean(ic()$datos$ls - ic()$datos$li)))

  output$ic_plot <- renderPlot({
    s <- ic()
    m <- m()
    k <- nrow(s$datos)
    g <- ggplot(s$datos, aes(y = muestra, colour = contiene)) +
      geom_segment(aes(x = li, xend = ls, yend = muestra), linewidth = if (k <= 100) 0.8 else 0.4) +
      geom_vline(xintercept = s$param, linetype = "dashed", colour = m$secundario, linewidth = 0.6) +
      scale_colour_manual(values = setNames(m$series[1:2], levels(s$datos$contiene)), name = NULL, drop = FALSE) +
      labs(
        x = "Intervalo de confianza", y = "Muestra",
        subtitle = paste("Línea discontinua: valor del parámetro", s$simbolo, "=", fmt(s$param, 2))
      ) +
      guides(colour = guide_legend(override.aes = list(linewidth = 1.5))) +
      tema(m) +
      theme(panel.grid.major.y = element_blank())
    if (k <= 100) g <- g + geom_point(aes(x = est), size = 1.5)
    g
  })

  output$ic_cob_plot <- renderPlot({
    s <- ic()
    m <- m()
    d <- data.frame(
      n = s$datos$muestra,
      cobertura = cumsum(s$datos$contiene == "Contiene el parámetro") / s$datos$muestra
    )
    ggplot(d, aes(n, cobertura)) +
      geom_hline(yintercept = s$conf, linetype = "dashed", colour = m$secundario, linewidth = 0.6) +
      geom_line(colour = m$series[1], linewidth = 0.8) +
      scale_y_continuous(labels = function(x) paste0(fmt(100 * x, 0), " %")) +
      labs(
        x = "Número de intervalos", y = "Proporción que contiene el parámetro",
        subtitle = "Línea discontinua: nivel de confianza"
      ) +
      tema(m)
  })

  # Contrastes de hipótesis ----------------------------------------------------

  ct_par <- reactive({
    s <- list(
      tipo = input$ct_tipo, n1 = input$ct_n1, n2 = input$ct_n2,
      alfa = input$ct_alfa, alt = input$ct_alt
    )
    if (s$tipo == "media") {
      req(is.finite(input$ct_mu1), is.finite(input$ct_mu2))
      c(s, list(
        mu1 = input$ct_mu1, mu2 = input$ct_mu2,
        sd1 = input$ct_sd1, sd2 = input$ct_sd2, var = input$ct_var
      ))
    } else {
      c(s, list(p1 = input$ct_p1, p2 = input$ct_p2))
    }
  })

  ct <- reactive({
    input$ct_sim
    s <- ct_par()
    list(
      s = s, d = diferencia(s),
      h0 = contraste(bajo_h0(s), input$ct_k),
      h1 = contraste(s, input$ct_k)
    )
  })

  ct_teo <- reactive({
    s <- ct_par()
    list(tipo1 = prob_rechazo(bajo_h0(s)), potencia = prob_rechazo(s))
  })

  output$ct_efecto_titulo <- renderText(
    if (ct_par()$tipo == "media") "Tamaño del efecto (d de Cohen)" else "Tamaño del efecto (h de Cohen)"
  )
  output$ct_efecto <- renderText(fmt(tamano_efecto(ct_par()), 2))
  output$ct_efecto_magnitud <- renderText(paste("Efecto", tolower(magnitud_efecto(tamano_efecto(ct_par())))))
  output$ct_tipo1 <- renderText(pct(mean(ct()$h0$rechaza)))
  output$ct_tipo1_teo <- renderText(
    paste0("α = ", pct(ct_par()$alfa, 0), " · Teórico: ", pct(ct_teo()$tipo1))
  )
  output$ct_tipo2 <- renderText(if (ct()$d == 0) "—" else pct(mean(!ct()$h1$rechaza)))
  output$ct_potencia <- renderText(if (ct()$d == 0) "—" else pct(mean(ct()$h1$rechaza)))
  output$ct_tipo2_teo <- renderText(
    if (diferencia(ct_par()) == 0) "H₁ no es cierta" else paste("Teórico:", pct(1 - ct_teo()$potencia))
  )
  output$ct_potencia_teo <- renderText(
    if (diferencia(ct_par()) == 0) "H₁ no es cierta" else paste("Teórica:", pct(ct_teo()$potencia))
  )

  output$ct_plot <- renderPlot({
    s <- ct_par()
    m <- m()
    modo <- input$ct_distr
    simbolo <- if (s$tipo == "media") "μ" else "p"
    niveles <- c(
      sprintf("Bajo H₀ (%s₁ = %s₂)", simbolo, simbolo),
      sprintf("Bajo H₁ (%s₁ − %s₂ = %s)", simbolo, simbolo, fmt(diferencia(s), 2))
    )
    colores <- setNames(m$series[1:2], niveles)
    crit <- criticos(s)
    teoricas <- list(distribucion(bajo_h0(s)), distribucion(s))
    # Se recortan las colas extremas (con muestras pequeñas el estadístico t
    # tiene colas muy pesadas); las densidades simuladas se calculan sobre el
    # total de contrastes.
    if (modo == "sim") {
      r <- ct()
      lim <- range(quantile(c(r$h0$t, r$h1$t), c(0.002, 0.998)), crit)
    } else {
      lim <- range(sapply(teoricas, function(f) f$q(c(0.002, 0.998))), crit)
    }
    zona <- switch(s$alt,
      bilateral = data.frame(xmin = c(-Inf, crit[2]), xmax = c(crit[1], Inf)),
      mayor = data.frame(xmin = crit, xmax = Inf),
      menor = data.frame(xmin = -Inf, xmax = crit)
    )
    g <- ggplot() +
      geom_rect(
        data = zona, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
        fill = m$series[8], alpha = 0.12
      )
    if (modo != "teo") {
      r <- ct()
      k <- length(r$h0$t)
      d <- data.frame(t = c(r$h0$t, r$h1$t), hip = factor(rep(niveles, each = k), levels = niveles))
      d <- d[d$t >= lim[1] & d$t <= lim[2], ]
      g <- g + geom_histogram(
        data = d, aes(t, y = after_stat(count / (k * width)), fill = hip),
        binwidth = diff(lim) / 60, boundary = 0, position = "identity",
        alpha = if (modo == "sim") 0.55 else 0.35
      )
    }
    if (modo != "sim") {
      x <- seq(lim[1], lim[2], length.out = 500)
      curvas <- data.frame(
        t = c(x, x), densidad = c(teoricas[[1]]$d(x), teoricas[[2]]$d(x)),
        hip = factor(rep(niveles, each = length(x)), levels = niveles)
      )
      if (modo == "teo") {
        g <- g + geom_area(
          data = curvas, aes(t, densidad, fill = hip),
          position = "identity", alpha = 0.3
        )
      }
      g <- g + geom_line(data = curvas, aes(t, densidad, colour = hip), linewidth = 1)
    }
    g +
      geom_vline(xintercept = crit, linetype = "dashed", colour = m$secundario, linewidth = 0.6) +
      scale_fill_manual(values = colores, name = NULL, aesthetics = c("fill", "colour")) +
      coord_cartesian(xlim = lim) +
      labs(
        x = "Estadístico de contraste", y = "Densidad",
        subtitle = paste0(
          "Zona sombreada: región de rechazo de H₀ (α = ", pct(s$alfa, 0), ").\n",
          "Área de H₀ en ella: error de tipo I; área de H₁: potencia."
        )
      ) +
      tema(m)
  })

  # Errores y potencia al variar el tamaño muestral o la variabilidad, con el
  # resto de parámetros fijos.
  ct_curva <- reactive({
    input$ct_sim
    s <- ct_par()
    k <- min(input$ct_k, 2000)
    if (input$ct_eje == "n") {
      x <- unique(round(exp(seq(log(2), log(500), length.out = 30))))
      r <- s$n2 / s$n1
      casos <- lapply(x, function(v) modifyList(s, list(n1 = v, n2 = max(2, round(r * v)))))
      etiqueta <- "Tamaño de la muestra 1 (n₁), con n₂/n₁ fijo"
      actual <- s$n1
    } else if (s$tipo == "media") {
      f <- seq(0.1, 3, length.out = 30)
      x <- f * s$sd1
      casos <- lapply(f, function(v) modifyList(s, list(sd1 = v * s$sd1, sd2 = v * s$sd2)))
      etiqueta <- "Desviación típica de la población 1 (σ₁), con σ₂/σ₁ fijo"
      actual <- s$sd1
    } else {
      # Se desplazan ambas proporciones manteniendo su diferencia: la varianza
      # p(1 − p) es máxima en p = 0,5 y disminuye hacia los extremos.
      dif <- s$p2 - s$p1
      x <- seq(max(0.01, 0.01 - dif), min(0.99, 0.99 - dif), length.out = 30)
      casos <- lapply(x, function(v) modifyList(s, list(p1 = v, p2 = min(max(v + dif, 0.01), 0.99))))
      etiqueta <- "Proporción de la población 1 (p₁), con p₂ − p₁ fija"
      actual <- s$p1
    }
    datos <- data.frame(
      x = x,
      tipo1 = vapply(casos, function(v) mean(contraste(bajo_h0(v), k)$rechaza), numeric(1)),
      potencia = vapply(casos, function(v) mean(contraste(v, k)$rechaza), numeric(1)),
      tipo1_teo = vapply(casos, function(v) prob_rechazo(bajo_h0(v)), numeric(1)),
      potencia_teo = vapply(casos, prob_rechazo, numeric(1))
    )
    list(datos = datos, etiqueta = etiqueta, actual = actual, s = s, log = input$ct_eje == "n")
  })

  output$ct_curva_titulo <- renderText(
    if (input$ct_eje == "n") {
      "Errores y potencia según el tamaño muestral"
    } else {
      "Errores y potencia según la variabilidad de las poblaciones"
    }
  )

  output$ct_curva_plot <- renderPlot({
    cv <- ct_curva()
    m <- m()
    series <- c("Error de tipo I", "Error de tipo II (β)", "Potencia (1 − β)")
    fuentes <- c("Simulado", "Teórico")
    d <- with(cv$datos, rbind(
      data.frame(x = x, valor = tipo1, serie = series[1], fuente = fuentes[1]),
      data.frame(x = x, valor = 1 - potencia, serie = series[2], fuente = fuentes[1]),
      data.frame(x = x, valor = potencia, serie = series[3], fuente = fuentes[1]),
      data.frame(x = x, valor = tipo1_teo, serie = series[1], fuente = fuentes[2]),
      data.frame(x = x, valor = 1 - potencia_teo, serie = series[2], fuente = fuentes[2]),
      data.frame(x = x, valor = potencia_teo, serie = series[3], fuente = fuentes[2])
    ))
    if (diferencia(cv$s) == 0) d <- d[d$serie == series[1], ]
    if (!input$ct_teorica) d <- d[d$fuente == fuentes[1], ]
    d$serie <- factor(d$serie, levels = series)
    d$fuente <- factor(d$fuente, levels = fuentes)
    g <- ggplot(d, aes(x, valor, colour = serie)) +
      geom_hline(yintercept = cv$s$alfa, linetype = "dashed", colour = m$secundario, linewidth = 0.6) +
      geom_vline(xintercept = cv$actual, linetype = "dotted", colour = m$secundario, linewidth = 0.6) +
      geom_line(aes(linetype = fuente, linewidth = fuente)) +
      geom_point(data = d[d$fuente == fuentes[1], ], size = 1.5) +
      scale_colour_manual(values = setNames(m$series[c(1, 3, 2)], series), name = NULL) +
      scale_linetype_manual(values = setNames(c("solid", "22"), fuentes), name = NULL) +
      scale_linewidth_manual(values = setNames(c(0.8, 0.6), fuentes), name = NULL) +
      scale_y_continuous(limits = c(0, 1), labels = function(x) pct(x, 0)) +
      labs(
        x = cv$etiqueta, y = "Probabilidad",
        subtitle = "Línea discontinua: α. Línea de puntos: valor actual."
      ) +
      guides(colour = guide_legend(order = 1)) +
      tema(m) +
      theme(legend.box = "vertical", legend.spacing.y = unit(0, "pt"))
    if (!input$ct_teorica) g <- g + guides(linetype = "none", linewidth = "none")
    if (cv$log) g <- g + scale_x_log10()
    g
  })
}

shinyApp(ui, server)
