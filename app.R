# Simulaciones de procesos estadísticos aleatorios
# - Ley de los grandes números
# - Aproximación de la distribución binomial mediante la Poisson
# - Teorema central del límite
# - Intervalos de confianza para la media y para proporciones
# - Contrastes de hipótesis sobre el parámetro de una población normal,
#   binomial o de Poisson
# - Comparación de dos poblaciones normales, binomiales o de Poisson

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

# Los parámetros de un contraste son una lista con el tipo ("media",
# "proporcion" o "poisson"), los tamaños muestrales n1 y n2, el nivel de
# significación alfa, la hipótesis alternativa alt ("bilateral", "mayor" o
# "menor") y los parámetros de cada población: mu1, mu2, sd1, sd2 y var
# ("welch" o "iguales") en las medias de poblaciones normales, p1 y p2 en las
# proporciones (n1 y n2 son entonces el número de ensayos) o l1 y l2 en las
# medias de poblaciones de Poisson.
parametros <- list(media = c("mu1", "mu2"), proporcion = c("p1", "p2"), poisson = c("l1", "l2"))

diferencia <- function(s) {
  k <- parametros[[s$tipo]]
  s[[k[1]]] - s[[k[2]]]
}

# Parámetros con θ₁ − θ₂ = d, cambiando el de la segunda población.
con_diferencia <- function(s, d) {
  k <- parametros[[s$tipo]]
  s[[k[2]]] <- s[[k[1]]] - d
  s
}

# Bajo la hipótesis nula la segunda población tiene la misma media (o
# proporción) que la primera; el resto de parámetros no cambia.
bajo_h0 <- function(s) con_diferencia(s, 0)

# Varianza de una observación en función del parámetro en las poblaciones
# binomiales (ensayos de Bernoulli) y de Poisson.
varianza <- function(s, theta) if (s$tipo == "proporcion") theta * (1 - theta) else theta

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
    # Totales de cada muestra: número de éxitos en las proporciones y suma de
    # las observaciones en las de Poisson. El estadístico usa la estimación
    # conjunta del parámetro bajo H0.
    if (s$tipo == "proporcion") {
      x1 <- rbinom(k, n1, s$p1)
      x2 <- rbinom(k, n2, s$p2)
    } else {
      x1 <- rpois(k, n1 * s$l1)
      x2 <- rpois(k, n2 * s$l2)
    }
    se <- sqrt(varianza(s, (x1 + x2) / (n1 + n2)) * (1 / n1 + 1 / n2))
    # Si todas las observaciones son iguales (todo éxitos o fracasos, o todo
    # ceros) no hay diferencia.
    t <- ifelse(se > 0, (x1 / n1 - x2 / n2) / se, 0)
    q <- qnorm
  }
  list(t = t, rechaza = rechazo(t, s$alt, s$alfa, q))
}

# Valores críticos del estadístico. En el contraste de Welch los grados de
# libertad dependen de cada muestra; se usan los de las varianzas poblacionales.
criticos <- function(s) {
  q <- if (s$tipo != "media") {
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
# - Proporciones y Poisson: aproximación normal de la diferencia de
#   proporciones (o medias) muestrales, con el error típico ponderado se0 en el
#   denominador del estadístico y el error típico real se1.
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
    k <- parametros[[s$tipo]]
    t1 <- s[[k[1]]]
    t2 <- s[[k[2]]]
    conjunto <- (s$n1 * t1 + s$n2 * t2) / (s$n1 + s$n2)
    se0 <- sqrt(varianza(s, conjunto) * (1 / s$n1 + 1 / s$n2))
    se1 <- sqrt(varianza(s, t1) / s$n1 + varianza(s, t2) / s$n2)
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
# si H0 es cierta y potencia si no lo es. Por defecto, para la comparación de
# dos poblaciones.
prob_rechazo <- function(s, dist = distribucion, crit = criticos) {
  f <- dist(s)$p
  crit <- crit(s)
  switch(s$alt,
    bilateral = f(crit[1]) + 1 - f(crit[2]),
    mayor = 1 - f(crit),
    menor = f(crit)
  )
}

# Contrastes sobre el parámetro de una población -------------------------------

# Los parámetros son una lista con el tipo ("media", "proporcion" o "poisson"),
# el número de observaciones obs, el tamaño muestral n (n = obs, salvo en las
# binomiales, donde es el número total de ensayos obs·m), el nivel de
# significación alfa, la hipótesis alternativa alt, el valor del parámetro en
# la población theta y el de la hipótesis nula theta0. En las medias de
# poblaciones normales, además, la desviación típica sd y sigma ("t" si es
# desconocida o "z" si es conocida).

# Simula k contrastes. En las normales se simulan directamente la media y la
# cuasivarianza muestrales; en las binomiales y de Poisson, el total de la
# muestra. El error típico de las binomiales y de Poisson usa el valor de H0.
contraste_una <- function(s, k) {
  n <- s$n
  if (s$tipo == "media") {
    m <- rnorm(k, s$theta, s$sd / sqrt(n))
    if (s$sigma == "z") {
      t <- (m - s$theta0) / (s$sd / sqrt(n))
    } else {
      v <- s$sd^2 * rchisq(k, n - 1) / (n - 1)
      t <- (m - s$theta0) / sqrt(v / n)
    }
  } else {
    x <- if (s$tipo == "proporcion") rbinom(k, n, s$theta) else rpois(k, n * s$theta)
    t <- (x / n - s$theta0) / sqrt(varianza(s, s$theta0) / n)
  }
  list(t = t, rechaza = rechazo(t, s$alt, s$alfa, cuantil_una(s)))
}

# Función cuantil del estadístico bajo H0: t con n − 1 grados de libertad en
# las medias con σ desconocida y normal estándar en los demás casos.
cuantil_una <- function(s) {
  if (s$tipo == "media" && s$sigma == "t") function(p) qt(p, s$n - 1) else qnorm
}

criticos_una <- function(s) {
  q <- cuantil_una(s)
  switch(s$alt,
    bilateral = c(-1, 1) * q(1 - s$alfa / 2),
    mayor = q(1 - s$alfa),
    menor = q(s$alfa)
  )
}

# Distribución teórica del estadístico (densidad d, distribución p y cuantil
# q). En las medias de poblaciones normales es exacta: normal con media
# (μ − μ₀)/(σ/√n) y desviación típica 1 si σ es conocida, o t no central con
# n − 1 grados de libertad y ese parámetro de no centralidad si no lo es. En
# las binomiales y de Poisson es la aproximación normal, con el error típico
# de H0 (se0) en el denominador del estadístico y el real (se1).
distribucion_una <- function(s) {
  if (s$tipo == "media") {
    ncp <- (s$theta - s$theta0) / (s$sd / sqrt(s$n))
    if (s$sigma == "z") {
      return(list(
        d = function(x) dnorm(x, ncp),
        p = function(x) pnorm(x, ncp),
        q = function(p) qnorm(p, ncp)
      ))
    }
    gl <- s$n - 1
    return(list(
      d = function(x) suppressWarnings(dt(x, gl, ncp)),
      p = function(x) suppressWarnings(pt(x, gl, ncp)),
      q = function(p) suppressWarnings(qt(p, gl, ncp))
    ))
  }
  se0 <- sqrt(varianza(s, s$theta0) / s$n)
  se1 <- sqrt(varianza(s, s$theta) / s$n)
  media <- (s$theta - s$theta0) / se0
  sd <- se1 / se0
  list(
    d = function(x) dnorm(x, media, sd),
    p = function(x) pnorm(x, media, sd),
    q = function(p) qnorm(p, media, sd)
  )
}

# Modelos de contraste --------------------------------------------------------

# Cada pestaña de contrastes usa un modelo con las funciones que dependen de
# si se compara una o dos poblaciones: simulación (sim), distribución teórica
# (dist), valores críticos (crit), parámetros bajo H0 (h0), diferencia con H0
# (dif), parámetros con una diferencia dada (con_dif), parámetros válidos
# (valido) y textos de los gráficos.

simbolos <- c(media = "μ", proporcion = "p", poisson = "λ")
num <- function(x) prettyNum(signif(x, 6), decimal.mark = ",")

modelo_una <- list(
  sim = contraste_una,
  dist = distribucion_una,
  crit = criticos_una,
  h0 = function(s) {
    s$theta <- s$theta0
    s
  },
  dif = function(s) s$theta - s$theta0,
  con_dif = function(s, d) {
    s$theta <- s$theta0 + d
    s
  },
  valido = function(s) {
    switch(s$tipo,
      media = TRUE,
      proporcion = s$theta > 0 && s$theta < 1,
      poisson = s$theta > 0
    )
  },
  etiquetas = function(s, h1) {
    sb <- simbolos[[s$tipo]]
    c(sprintf("%s = %s", sb, num(s$theta0)), sprintf("%s = %s", sb, num(h1$theta)))
  },
  subtitulo = function(s) {
    sb <- simbolos[[s$tipo]]
    sprintf("Población: %s = %s (H₀: %s = %s)", sb, num(s$theta), sb, num(s$theta0))
  },
  # Fórmula del estadístico en LaTeX, con su distribución bajo H0.
  latex = function(s) {
    # En las binomiales el tamaño es el número total de ensayos n·m.
    n <- if (s$m > 1) r"(n \cdot m)" else "n"
    switch(s$tipo,
      media = if (s$sigma == "t") {
        r"(T = \frac{\bar{X} - \mu_0}{S/\sqrt{n}} \sim t_{n-1} \quad \text{bajo } H_0)"
      } else {
        r"(Z = \frac{\bar{X} - \mu_0}{\sigma/\sqrt{n}} \sim N(0, 1) \quad \text{bajo } H_0)"
      },
      proporcion = sprintf(
        r"(Z = \frac{\hat{P} - p_0}{\sqrt{\dfrac{p_0(1 - p_0)}{%s}}} \sim N(0, 1) \quad \text{aprox. bajo } H_0)", n
      ),
      poisson = r"(Z = \frac{\bar{X} - \lambda_0}{\sqrt{\dfrac{\lambda_0}{n}}} \sim N(0, 1) \quad \text{aprox. bajo } H_0)"
    )
  },
  fuera_de_rango = "Con esta diferencia mínima el parámetro bajo H₁ queda fuera de su rango."
)

modelo_dos <- list(
  sim = contraste,
  dist = distribucion,
  crit = criticos,
  h0 = bajo_h0,
  dif = diferencia,
  con_dif = con_diferencia,
  valido = function(s) {
    switch(s$tipo,
      media = TRUE,
      proporcion = s$p2 > 0 && s$p2 < 1,
      poisson = s$l2 > 0
    )
  },
  etiquetas = function(s, h1) {
    sb <- simbolos[[s$tipo]]
    c(sprintf("%s₁ = %s₂", sb, sb), sprintf("%s₁ − %s₂ = %s", sb, sb, fmt(diferencia(h1), 2)))
  },
  subtitulo = function(s) {
    sb <- simbolos[[s$tipo]]
    sprintf("Poblaciones: %s₁ − %s₂ = %s", sb, sb, fmt(diferencia(s), 2))
  },
  # Fórmula del estadístico en LaTeX, con su distribución bajo H0.
  latex = function(s) {
    # En las binomiales los tamaños son los números totales de ensayos nᵢ·m.
    n <- if (isTRUE(s$m > 1)) c(r"(n_1 \cdot m)", r"(n_2 \cdot m)") else c("n_1", "n_2")
    switch(s$tipo,
      media = if (s$var == "iguales") {
        paste(
          r"(T = \frac{\bar{X}_1 - \bar{X}_2}{S_p\sqrt{\dfrac{1}{n_1} + \dfrac{1}{n_2}}})",
          r"(\sim t_{n_1 + n_2 - 2} \quad \text{bajo } H_0)"
        )
      } else {
        paste(
          r"(T = \frac{\bar{X}_1 - \bar{X}_2}{\sqrt{\dfrac{S_1^2}{n_1} + \dfrac{S_2^2}{n_2}}})",
          r"(\sim t_{\nu} \quad \text{aprox. bajo } H_0)"
        )
      },
      proporcion = sprintf(
        paste(
          r"(Z = \frac{\hat{P}_1 - \hat{P}_2}{\sqrt{\hat{P}(1 - \hat{P})\left(\dfrac{1}{%s} + \dfrac{1}{%s}\right)}})",
          r"(\sim N(0, 1) \quad \text{aprox. bajo } H_0)"
        ),
        n[1], n[2]
      ),
      poisson = paste(
        r"(Z = \frac{\bar{X}_1 - \bar{X}_2}{\sqrt{\hat{\lambda}\left(\dfrac{1}{n_1} + \dfrac{1}{n_2}\right)}})",
        r"(\sim N(0, 1) \quad \text{aprox. bajo } H_0)"
      )
    )
  },
  fuera_de_rango = "Con esta diferencia mínima el parámetro de la población 2 bajo H₁ queda fuera de su rango."
)

prob_modelo <- function(s, modelo) prob_rechazo(s, modelo$dist, modelo$crit)

# Parámetros bajo H1 con la diferencia mínima δ (s$delta) en el sentido de la
# alternativa. En la bilateral se toma el sentido de la diferencia real, salvo
# que deje el parámetro fuera de su rango. Devuelve NULL si no hay ninguno
# válido.
h1_minima <- function(s, modelo) {
  signo <- switch(s$alt, mayor = 1, menor = -1, bilateral = if (modelo$dif(s) < 0) -1 else 1)
  h1 <- modelo$con_dif(s, signo * s$delta)
  if (s$alt == "bilateral" && !modelo$valido(h1)) h1 <- modelo$con_dif(s, -signo * s$delta)
  if (modelo$valido(h1)) h1
}

# Rango de los tamaños muestrales según la distribución de las poblaciones.
rangos_n <- list(normal = c(5, 100), binomial = c(10, 1000), poisson = c(5, 100))

# Opacidad de las barras del histograma del estadístico, semitransparentes
# para que se vean las áreas de los riesgos α y β, y de esas áreas.
opacidad_barras <- 0.45
opacidad_areas <- 0.3

# Tema de un recuadro de valores con el mismo color que su elemento del
# gráfico: el color con su opacidad sobre el fondo. El texto es el color de
# texto o el de fondo de la app, el que más contraste tenga.
tema_caja <- function(m, color, opacidad) {
  bg <- opacidad * col2rgb(color) + (1 - opacidad) * col2rgb(m$superficie)
  luminancia <- function(rgb) {
    c <- rgb / 255
    c <- ifelse(c <= 0.03928, c / 12.92, ((c + 0.055) / 1.055)^2.4)
    sum(c(0.2126, 0.7152, 0.0722) * c)
  }
  contraste <- function(fg) {
    l <- sort(c(luminancia(bg), luminancia(col2rgb(fg))), decreasing = TRUE)
    (l[1] + 0.05) / (l[2] + 0.05)
  }
  fg <- c(m$texto, m$superficie)
  value_box_theme(
    bg = rgb(t(bg), maxColorValue = 255),
    fg = fg[which.max(vapply(fg, contraste, numeric(1)))]
  )
}

# Servidor común de las pestañas de contrastes, con los controles y salidas
# con el prefijo dado: deslizadores de tamaño muestral según la distribución,
# simulación al pulsar el botón, recuadros de aceptaciones, rechazos y riesgos
# α y β, y gráfico de la distribución del estadístico. par es la reacción con
# los parámetros del contraste (incluida la diferencia mínima delta) y ids_n
# los deslizadores de los tamaños muestrales. Devuelve las reacciones con los
# parámetros, los parámetros bajo H1, la clave de la simulación (parámetros y
# número de simulaciones) y la simulación (NULL si no se ha simulado con los
# parámetros actuales).
servidor_contraste <- function(input, output, session, prefijo, par, modelo, ids_n, m) {
  id <- function(x) paste0(prefijo, "_", x)
  entrada <- function(x) input[[id(x)]]

  # Al cambiar de distribución se ajustan los deslizadores de los tamaños
  # muestrales, conservando su valor si cabe.
  observeEvent(entrada("dist"), {
    r <- rangos_n[[entrada("dist")]]
    for (i in ids_n) {
      updateSliderInput(session, id(i), min = r[1], max = r[2], value = min(max(entrada(i), r[1]), r[2]))
    }
  }, ignoreInit = TRUE)

  h1 <- reactive({
    s <- par()
    validate(need(isTRUE(s$delta > 0), "La diferencia mínima a detectar debe ser positiva."))
    h <- h1_minima(s, modelo)
    validate(need(!is.null(h), modelo$fuera_de_rango))
    h
  })

  # Simulación al pulsar el botón. Se guarda con los parámetros usados y solo
  # se muestra mientras coincidan con los actuales: al cambiar cualquier
  # parámetro desaparece hasta que se vuelva a simular.
  clave <- reactive(c(par(), k = entrada("k")))
  simulado <- reactiveVal(NULL)
  observeEvent(entrada("sim"), {
    req(entrada("sim") > 0)
    s <- clave()
    simulado(list(clave = s, r = modelo$sim(s, s$k)))
  })
  sim <- reactive({
    x <- simulado()
    if (!is.null(x) && identical(x$clave, clave())) x$r
  })

  # ¿Es cierta H0 en las poblaciones simuladas? En las alternativas
  # unilaterales H0 incluye también las diferencias en el sentido contrario.
  h0_cierta <- reactive({
    d <- modelo$dif(par())
    switch(par()$alt, bilateral = d == 0, mayor = d <= 0, menor = d >= 0)
  })

  output[[id("valores")]] <- renderUI({
    m <- m()
    fila_valores(
      value_box(
        "Aceptaciones de H₀", textOutput(id("acepta")),
        p(textOutput(id("acepta_nota"), inline = TRUE)),
        theme = tema_caja(m, m$series[3], opacidad_barras)
      ),
      value_box(
        "Rechazos de H₀", textOutput(id("rechaza")),
        p(textOutput(id("rechaza_nota"), inline = TRUE)),
        theme = tema_caja(m, m$series[8], opacidad_barras)
      ),
      value_box(
        "Riesgo α", textOutput(id("alfa")),
        p(textOutput(id("alfa_nota"), inline = TRUE)),
        theme = tema_caja(m, m$series[1], opacidad_areas)
      ),
      value_box(
        "Riesgo β (con la diferencia mínima)", textOutput(id("beta")),
        p(textOutput(id("potencia"), inline = TRUE)),
        theme = tema_caja(m, m$series[2], opacidad_areas)
      )
    )
  })

  # Número y porcentaje de rechazos (o de aceptaciones) de la simulación.
  frecuencia <- function(rechazos) {
    r <- sim()
    if (is.null(r)) return("—")
    x <- r$rechaza == rechazos
    sprintf("%d (%s)", sum(x), pct(mean(x)))
  }
  nota <- function(cierta, falsa) {
    if (is.null(sim())) "Pulsa «Simular»" else if (h0_cierta()) cierta else falsa
  }
  output[[id("acepta")]] <- renderText(frecuencia(FALSE))
  output[[id("rechaza")]] <- renderText(frecuencia(TRUE))
  output[[id("acepta_nota")]] <- renderText(
    nota("Decisiones correctas (H₀ cierta)", "Errores de tipo II (H₀ falsa)")
  )
  output[[id("rechaza_nota")]] <- renderText(
    nota("Errores de tipo I (H₀ cierta)", "Decisiones correctas (H₀ falsa)")
  )
  output[[id("alfa")]] <- renderText(pct(prob_modelo(modelo$h0(par()), modelo)))
  output[[id("alfa_nota")]] <- renderText(paste("Nivel de significación:", pct(par()$alfa, 0)))
  output[[id("beta")]] <- renderText(pct(1 - prob_modelo(h1(), modelo)))
  output[[id("potencia")]] <- renderText(paste("Potencia (1 − β):", pct(prob_modelo(h1(), modelo))))

  # Fórmula del estadístico, que MathJax compone cada vez que cambia.
  output[[id("formula")]] <- renderUI(
    withMathJax(div(class = "formula-contraste", paste0("$$", modelo$latex(par()), "$$")))
  )

  # Distribuciones teóricas del estadístico bajo H0 y bajo H1.
  teoricas <- reactive(list(modelo$dist(modelo$h0(par())), modelo$dist(h1())))

  # Rango del eje X del estadístico: el de las distribuciones teóricas, los
  # valores críticos y, si se ha simulado, los valores simulados. Se recortan
  # las colas extremas; las frecuencias se calculan sobre el total de
  # simulaciones.
  lim_auto <- reactive({
    lim <- range(sapply(teoricas(), function(f) f$q(c(0.002, 0.998))), modelo$crit(par()))
    if (!is.null(sim())) lim <- range(lim, quantile(sim()$t, c(0.002, 0.998)))
    lim
  })

  # Al fijar la escala se guarda el rango del momento; al soltarla se vuelve
  # al rango automático.
  lim_fijo <- reactiveVal(NULL)
  observeEvent(entrada("eje_fijo"), {
    lim_fijo(if (entrada("eje_fijo")) tryCatch(lim_auto(), error = function(e) NULL))
  })

  output[[id("plot")]] <- renderPlot({
    s <- par()
    h1 <- h1()
    teoricas <- teoricas()
    r <- sim()
    m <- m()
    crit <- modelo$crit(s)
    lim <- lim_fijo()
    if (is.null(lim)) lim <- lim_auto()
    g <- ggplot()
    if (!is.null(r)) {
      k <- length(r$t)
      decisiones <- c("Aceptación de H₀", "Rechazo de H₀")
      d <- data.frame(
        t = r$t,
        decision = factor(ifelse(r$rechaza, decisiones[2], decisiones[1]), levels = decisiones)
      )
      d <- d[d$t >= lim[1] & d$t <= lim[2], ]
      # Ancho de clase con un número entero de clases entre los dos valores
      # críticos, para que ninguna clase quede a ambos lados de uno de ellos.
      ancho <- diff(lim) / 60
      if (length(crit) == 2) ancho <- diff(crit) / max(1, round(diff(crit) / ancho))
      g <- g +
        geom_histogram(
          data = d, aes(t, y = after_stat(count / (k * width)), fill = decision),
          binwidth = ancho, boundary = crit[1], colour = m$superficie, linewidth = 0.2,
          alpha = opacidad_barras,
          # Muestra en la leyenda también la decisión que no se haya dado nunca.
          show.legend = c(fill = TRUE, colour = FALSE)
        ) +
        scale_fill_manual(values = setNames(m$series[c(3, 8)], decisiones), name = NULL, drop = FALSE)
    }
    x <- seq(lim[1], lim[2], length.out = 500)
    etiquetas <- modelo$etiquetas(s, h1)
    curvas <- c(
      sprintf("Teórica bajo H₀ (%s), área sombreada: α", etiquetas[1]),
      sprintf("Teórica bajo H₁ (%s), área sombreada: β", etiquetas[2])
    )
    dc <- data.frame(
      t = c(x, x), densidad = c(teoricas[[1]]$d(x), teoricas[[2]]$d(x)),
      curva = factor(rep(curvas, each = length(x)), levels = curvas)
    )
    # Riesgo α: área bajo la densidad de H0 en la región de rechazo. Riesgo β:
    # área bajo la densidad de H1 en la región de aceptación. Cada intervalo
    # se recorta a los límites del gráfico.
    rechazo <- switch(s$alt,
      bilateral = list(c(-Inf, crit[1]), c(crit[2], Inf)),
      mayor = list(c(crit, Inf)),
      menor = list(c(-Inf, crit))
    )
    aceptacion <- switch(s$alt,
      bilateral = list(crit),
      mayor = list(c(-Inf, crit)),
      menor = list(c(crit, Inf))
    )
    areas <- function(f, intervalos, nombre) {
      do.call(rbind, lapply(seq_along(intervalos), function(i) {
        a <- max(intervalos[[i]][1], lim[1])
        b <- min(intervalos[[i]][2], lim[2])
        if (a >= b) return(NULL)
        x <- seq(a, b, length.out = 200)
        data.frame(t = x, densidad = f(x), trozo = paste(nombre, i))
      }))
    }
    da <- areas(teoricas[[1]]$d, rechazo, "alfa")
    db <- areas(teoricas[[2]]$d, aceptacion, "beta")
    if (!is.null(da)) {
      g <- g + geom_ribbon(
        data = da, aes(t, ymin = 0, ymax = densidad, group = trozo),
        fill = m$series[1], alpha = opacidad_areas
      )
    }
    if (!is.null(db)) {
      g <- g + geom_ribbon(
        data = db, aes(t, ymin = 0, ymax = densidad, group = trozo),
        fill = m$series[2], alpha = opacidad_areas
      )
    }
    g +
      geom_line(data = dc, aes(t, densidad, colour = curva), linewidth = 1) +
      geom_vline(xintercept = crit, linetype = "dashed", colour = m$secundario, linewidth = 0.6) +
      scale_colour_manual(values = setNames(m$series[1:2], curvas), name = NULL) +
      guides(fill = guide_legend(order = 1), colour = guide_legend(ncol = 1, order = 2)) +
      # El eje Y se ajusta a las curvas teóricas: con pocas simulaciones las
      # barras del histograma son muy altas y las aplastarían, así que se
      # cortan por arriba.
      coord_cartesian(xlim = lim, ylim = c(0, 1.5 * max(dc$densidad))) +
      labs(
        x = "Estadístico de contraste", y = "Densidad",
        subtitle = paste0(
          modelo$subtitulo(s), ". Línea discontinua: valores críticos.",
          if (is.null(r)) "\nPulsa «Simular» para ver la distribución de frecuencias simulada."
        )
      ) +
      tema(m) +
      theme(legend.box = "vertical", legend.spacing.y = unit(0, "pt"))
  })

  list(par = par, h1 = h1, clave = clave, sim = sim)
}

# Interfaz ---------------------------------------------------------------------

boton_simular <- function(id, etiqueta = "Simular de nuevo") {
  actionButton(id, etiqueta, icon = icon("dice"), class = "btn-primary")
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
/* Fórmula del estadístico de contraste; se desplaza en horizontal si no cabe. */
.formula-contraste { overflow-x: auto; overflow-y: hidden; }
.formula-contraste mjx-container[display='true'], .formula-contraste .MathJax_Display {
  margin: 0.25rem 0 !important;
}
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
  "Contrastes de 1 población",
  layout_sidebar(
    sidebar = sidebar(
      width = 300,
      title = "Parámetros",
      selectInput("ct_dist", "Distribución de la población",
        choices = c("Normal" = "normal", "Binomial" = "binomial", "Poisson" = "poisson")
      ),
      conditionalPanel(
        "input.ct_dist == 'normal'",
        numericInput("ct_mu0", "Media bajo H₀ (μ₀)", value = 100),
        numericInput("ct_mu", "Media de la población (μ)", value = 105),
        numericInput("ct_sd", "Desviación típica de la población (σ)", value = 15, min = 0.01),
        radioButtons("ct_sigma", "Desviación típica poblacional",
          choices = c("Desconocida (contraste t)" = "t", "Conocida (contraste z)" = "z")
        ),
        numericInput("ct_delta_norm", "Diferencia mínima a detectar (δ = |μ − μ₀|)", value = 5, min = 0)
      ),
      conditionalPanel(
        "input.ct_dist == 'binomial'",
        sliderInput("ct_m", "Número de ensayos de cada observación (m)", min = 1, max = 50, value = 1, step = 1),
        sliderInput("ct_p0", "Probabilidad de éxito bajo H₀ (p₀)", min = 0.01, max = 0.99, value = 0.3, step = 0.01),
        sliderInput("ct_p", "Probabilidad de éxito en la población (p)", min = 0.01, max = 0.99, value = 0.4, step = 0.01),
        sliderInput("ct_delta_binom", "Diferencia mínima a detectar (δ = |p − p₀|)",
          min = 0.01, max = 0.5, value = 0.1, step = 0.01
        )
      ),
      conditionalPanel(
        "input.ct_dist == 'poisson'",
        sliderInput("ct_l0", "Media bajo H₀ (λ₀)", min = 0.1, max = 20, value = 4, step = 0.1),
        sliderInput("ct_l", "Media de la población (λ)", min = 0.1, max = 20, value = 5, step = 0.1),
        sliderInput("ct_delta_pois", "Diferencia mínima a detectar (δ = |λ − λ₀|)",
          min = 0.1, max = 10, value = 1, step = 0.1
        )
      ),
      # El rango del tamaño muestral depende de la distribución (ver rangos_n);
      # este es el de la normal.
      sliderInput("ct_n", "Tamaño muestral (n)", min = 5, max = 100, value = 30, step = 1),
      selectInput("ct_alt", "Hipótesis alternativa (H₁)",
        choices = c(
          "Bilateral: parámetro ≠ valor de H₀" = "bilateral",
          "Unilateral: parámetro > valor de H₀" = "mayor",
          "Unilateral: parámetro < valor de H₀" = "menor"
        )
      ),
      sliderInput("ct_alfa", "Nivel de significación (α)", min = 0.01, max = 0.20, value = 0.05, step = 0.01),
      sliderInput("ct_k", "Número de simulaciones", min = 1, max = 10000, value = 1000, step = 1),
      radioButtons("ct_eje", "En el gráfico de errores, variar",
        choices = c("El tamaño muestral" = "n", "La variabilidad de la población" = "var")
      ),
      checkboxInput("ct_eje_fijo", "Fijar la escala del eje X en el rango actual", value = FALSE),
      boton_simular("ct_sim", "Simular"),
      helpText(
        "En cada simulación se extrae una muestra de la población y se realiza",
        "el contraste. La simulación se borra al cambiar cualquier parámetro y hay",
        "que volver a pulsar «Simular». Las distribuciones teóricas se calculan",
        "bajo H₀ y bajo H₁ con la diferencia mínima δ en el sentido de la",
        "alternativa (en la bilateral, en el de la diferencia real con H₀)."
      ),
      helpText(
        "Contrastes: t (o z si σ es conocida) para la media de una población",
        "normal y z con el error típico de H₀ para las binomiales y las de Poisson."
      )
    ),
    uiOutput("ct_valores"),
    layout_columns(
      col_widths = c(6, 6),
      card(
        full_screen = TRUE,
        card_header("Distribución del estadístico de contraste"),
        uiOutput("ct_formula"),
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

ui_cd <- nav_panel(
  "Contrastes de 2 poblaciones",
  layout_sidebar(
    sidebar = sidebar(
      width = 300,
      title = "Parámetros",
      selectInput("cd_dist", "Distribución de las poblaciones",
        choices = c("Normal" = "normal", "Binomial" = "binomial", "Poisson" = "poisson")
      ),
      conditionalPanel(
        "input.cd_dist == 'normal'",
        numericInput("cd_mu1", "Media de la población 1 (μ₁)", value = 100),
        numericInput("cd_sd1", "Desviación típica de la población 1 (σ₁)", value = 15, min = 0.01),
        numericInput("cd_mu2", "Media de la población 2 (μ₂)", value = 105),
        numericInput("cd_sd2", "Desviación típica de la población 2 (σ₂)", value = 15, min = 0.01),
        radioButtons("cd_var", "Varianzas poblacionales (desconocidas)",
          choices = c(
            "Distintas (t de Welch)" = "welch",
            "Iguales (t de Student con varianza ponderada)" = "iguales"
          )
        ),
        numericInput("cd_delta_norm", "Diferencia mínima a detectar (δ = |μ₁ − μ₂|)", value = 5, min = 0)
      ),
      conditionalPanel(
        "input.cd_dist == 'binomial'",
        sliderInput("cd_m", "Número de ensayos de cada observación (m)", min = 1, max = 50, value = 1, step = 1),
        sliderInput("cd_p1", "Probabilidad de éxito en la población 1 (p₁)", min = 0.01, max = 0.99, value = 0.3, step = 0.01),
        sliderInput("cd_p2", "Probabilidad de éxito en la población 2 (p₂)", min = 0.01, max = 0.99, value = 0.4, step = 0.01),
        sliderInput("cd_delta_binom", "Diferencia mínima a detectar (δ = |p₁ − p₂|)",
          min = 0.01, max = 0.5, value = 0.1, step = 0.01
        )
      ),
      conditionalPanel(
        "input.cd_dist == 'poisson'",
        sliderInput("cd_l1", "Media de la población 1 (λ₁)", min = 0.1, max = 20, value = 4, step = 0.1),
        sliderInput("cd_l2", "Media de la población 2 (λ₂)", min = 0.1, max = 20, value = 5, step = 0.1),
        sliderInput("cd_delta_pois", "Diferencia mínima a detectar (δ = |λ₁ − λ₂|)",
          min = 0.1, max = 10, value = 1, step = 0.1
        )
      ),
      # El rango de los tamaños muestrales depende de la distribución (ver
      # rangos_n en el servidor); estos son los de la normal.
      sliderInput("cd_n1", "Tamaño de la muestra 1 (n₁)", min = 5, max = 100, value = 50, step = 1),
      sliderInput("cd_n2", "Tamaño de la muestra 2 (n₂)", min = 5, max = 100, value = 50, step = 1),
      selectInput("cd_alt", "Hipótesis alternativa (H₁)",
        choices = c(
          "Bilateral: población 1 ≠ población 2" = "bilateral",
          "Unilateral: población 1 > población 2" = "mayor",
          "Unilateral: población 1 < población 2" = "menor"
        )
      ),
      sliderInput("cd_alfa", "Nivel de significación (α)", min = 0.01, max = 0.20, value = 0.05, step = 0.01),
      sliderInput("cd_k", "Número de simulaciones", min = 1, max = 10000, value = 1000, step = 1),
      checkboxInput("cd_eje_fijo", "Fijar la escala del eje X en el rango actual", value = FALSE),
      boton_simular("cd_sim", "Simular"),
      helpText(
        "En cada simulación se extrae una muestra de cada población y se realiza",
        "el contraste. La simulación se borra al cambiar cualquier parámetro y hay",
        "que volver a pulsar «Simular». Las distribuciones teóricas se calculan bajo H₀ (con el",
        "parámetro de la población 2 igual al de la 1) y bajo H₁ con la diferencia",
        "mínima δ en el sentido de la alternativa (en la bilateral, en el de la",
        "diferencia real entre las poblaciones)."
      ),
      helpText(
        "Contrastes: t para las medias de poblaciones normales y z con la",
        "estimación conjunta del parámetro para las binomiales y las de Poisson."
      )
    ),
    uiOutput("cd_valores"),
    layout_columns(
      col_widths = c(8, 4),
      card(
        full_screen = TRUE,
        card_header("Distribución de frecuencias del estadístico de contraste"),
        uiOutput("cd_formula"),
        plotOutput("cd_plot", height = "450px")
      ),
      card(
        full_screen = TRUE,
        card_header("Distribución de las poblaciones"),
        plotOutput("cd_pob_plot", height = "450px")
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
  ui_cd,
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
    s <- list(obs = input$ct_n, n = input$ct_n, m = 1, alfa = input$ct_alfa, alt = input$ct_alt)
    switch(input$ct_dist,
      normal = {
        req(is.finite(input$ct_mu0), is.finite(input$ct_mu))
        validate(need(isTRUE(input$ct_sd > 0), "La desviación típica debe ser positiva."))
        c(s, list(
          tipo = "media", theta = input$ct_mu, theta0 = input$ct_mu0,
          sd = input$ct_sd, sigma = input$ct_sigma, delta = input$ct_delta_norm
        ))
      },
      binomial = {
        # Cada observación son m ensayos de Bernoulli, así que se contrasta la
        # proporción de éxitos en n·m ensayos.
        s$m <- input$ct_m
        s$n <- s$obs * s$m
        c(s, list(tipo = "proporcion", theta = input$ct_p, theta0 = input$ct_p0, delta = input$ct_delta_binom))
      },
      poisson = c(s, list(tipo = "poisson", theta = input$ct_l, theta0 = input$ct_l0, delta = input$ct_delta_pois))
    )
  })

  ct <- servidor_contraste(input, output, session, "ct", ct_par, modelo_una, "n", m)

  # Casos del gráfico de errores y potencia: cambia el tamaño muestral o la
  # variabilidad y el resto de parámetros se mantiene.
  casos_curva <- function(s, eje, dist) {
    if (eje == "n") {
      r <- rangos_n[[dist]]
      x <- unique(round(exp(seq(log(r[1]), log(r[2]), length.out = 30))))
      list(
        x = x, casos = lapply(x, function(v) modifyList(s, list(obs = v, n = v * s$m))),
        etiqueta = "Tamaño muestral (n)", actual = s$obs, log = TRUE
      )
    } else if (s$tipo == "media") {
      f <- seq(0.1, 3, length.out = 30)
      list(
        x = f * s$sd, casos = lapply(f, function(v) modifyList(s, list(sd = v * s$sd))),
        etiqueta = "Desviación típica de la población (σ)", actual = s$sd, log = FALSE
      )
    } else {
      # Se desplazan el valor de H0 y el de la población manteniendo su
      # diferencia: la varianza de una observación, p(1 − p) o λ, cambia con
      # el parámetro.
      dif <- s$theta - s$theta0
      if (s$tipo == "proporcion") {
        x <- seq(0.01, 0.99, length.out = 30)
        etiqueta <- "Proporción bajo H₀ (p₀), con p − p₀ fija"
      } else {
        x <- seq(0.1, max(20, 2 * s$theta0), length.out = 30)
        etiqueta <- "Media bajo H₀ (λ₀), con λ − λ₀ fija"
      }
      list(
        x = x, casos = lapply(x, function(v) modifyList(s, list(theta0 = v, theta = v + dif))),
        etiqueta = etiqueta, actual = s$theta0, log = FALSE
      )
    }
  }

  # Errores de tipo I (bajo H0) y potencia (bajo H1 con la diferencia mínima)
  # de cada caso, teóricos o simulados con k contrastes. Donde la diferencia
  # mínima deja el parámetro fuera de su rango la potencia es NA.
  errores_curva <- function(cs, k = NULL) {
    prob <- if (is.null(k)) {
      function(s) prob_modelo(s, modelo_una)
    } else {
      function(s) mean(contraste_una(s, k)$rechaza)
    }
    data.frame(
      x = cs$x,
      tipo1 = vapply(cs$casos, function(v) prob(modelo_una$h0(v)), numeric(1)),
      potencia = vapply(cs$casos, function(v) {
        h1 <- h1_minima(v, modelo_una)
        if (is.null(h1)) NA_real_ else prob(h1)
      }, numeric(1))
    )
  }

  ct_curva <- reactive({
    ct$h1()
    cs <- casos_curva(ct_par(), input$ct_eje, input$ct_dist)
    c(cs, list(datos = errores_curva(cs)))
  })

  # Al simular se calculan también los errores simulados del gráfico, para los
  # dos ejes posibles (con un máximo de 2000 contrastes por punto). Como la
  # simulación principal, solo se muestran con los parámetros con que se
  # hicieron.
  ct_curva_sim <- reactiveVal(NULL)
  observeEvent(input$ct_sim, {
    req(input$ct_sim > 0)
    s <- ct_par()
    k <- min(input$ct_k, 2000)
    ct_curva_sim(list(
      clave = ct$clave(),
      datos = lapply(c(n = "n", var = "var"), function(eje) errores_curva(casos_curva(s, eje, input$ct_dist), k))
    ))
  })

  output$ct_curva_titulo <- renderText(
    if (input$ct_eje == "n") {
      "Errores y potencia según el tamaño muestral"
    } else {
      "Errores y potencia según la variabilidad de la población"
    }
  )

  output$ct_curva_plot <- renderPlot({
    cv <- ct_curva()
    m <- m()
    simulado <- ct_curva_sim()
    simulado <- if (!is.null(simulado) && identical(simulado$clave, ct$clave())) simulado$datos[[input$ct_eje]]
    series <- c("Error de tipo I", "Error de tipo II (β)", "Potencia (1 − β)")
    fuentes <- c("Simulado", "Teórico")
    largo <- function(datos, fuente) {
      with(datos, rbind(
        data.frame(x = x, valor = tipo1, serie = series[1], fuente = fuente),
        data.frame(x = x, valor = 1 - potencia, serie = series[2], fuente = fuente),
        data.frame(x = x, valor = potencia, serie = series[3], fuente = fuente)
      ))
    }
    d <- largo(cv$datos, fuentes[2])
    if (!is.null(simulado)) d <- rbind(largo(simulado, fuentes[1]), d)
    d$serie <- factor(d$serie, levels = series)
    d$fuente <- factor(d$fuente, levels = fuentes)
    g <- ggplot(d, aes(x, valor, colour = serie)) +
      geom_hline(yintercept = ct_par()$alfa, linetype = "dashed", colour = m$secundario, linewidth = 0.6) +
      geom_vline(xintercept = cv$actual, linetype = "dotted", colour = m$secundario, linewidth = 0.6) +
      geom_line(aes(linetype = fuente, linewidth = fuente), na.rm = TRUE) +
      scale_colour_manual(values = setNames(m$series[c(1, 3, 2)], series), name = NULL) +
      # Las curvas teóricas son discontinuas si se superponen a las simuladas.
      scale_linetype_manual(
        values = setNames(c("solid", if (is.null(simulado)) "solid" else "22"), fuentes),
        name = NULL, drop = FALSE
      ) +
      scale_linewidth_manual(values = setNames(c(0.8, 0.6), fuentes), name = NULL, drop = FALSE) +
      scale_y_continuous(limits = c(0, 1), labels = function(x) pct(x, 0)) +
      labs(
        x = cv$etiqueta, y = "Probabilidad",
        subtitle = paste0(
          "Error de tipo II y potencia con la diferencia mínima δ.\n",
          "Línea discontinua: α. Línea de puntos: valor actual.",
          if (is.null(simulado)) "\nPulsa «Simular» para añadir los valores simulados."
        )
      ) +
      guides(colour = guide_legend(order = 1)) +
      tema(m) +
      theme(legend.box = "vertical", legend.spacing.y = unit(0, "pt"))
    if (!is.null(simulado)) {
      g <- g + geom_point(data = d[d$fuente == fuentes[1], ], size = 1.5, na.rm = TRUE)
    } else {
      g <- g + guides(linetype = "none", linewidth = "none")
    }
    if (cv$log) g <- g + scale_x_log10()
    g
  })

  # Comparación de dos poblaciones ---------------------------------------------

  cd_par <- reactive({
    s <- list(n1 = input$cd_n1, n2 = input$cd_n2, alfa = input$cd_alfa, alt = input$cd_alt)
    switch(input$cd_dist,
      normal = {
        req(is.finite(input$cd_mu1), is.finite(input$cd_mu2))
        validate(need(
          isTRUE(input$cd_sd1 > 0 && input$cd_sd2 > 0),
          "Las desviaciones típicas deben ser positivas."
        ))
        c(s, list(
          tipo = "media", mu1 = input$cd_mu1, mu2 = input$cd_mu2,
          sd1 = input$cd_sd1, sd2 = input$cd_sd2, var = input$cd_var,
          delta = input$cd_delta_norm
        ))
      },
      binomial = {
        # Cada observación son m ensayos de Bernoulli, así que se comparan las
        # proporciones de éxitos en n·m ensayos.
        s$n1 <- s$n1 * input$cd_m
        s$n2 <- s$n2 * input$cd_m
        c(s, list(
          tipo = "proporcion", m = input$cd_m, p1 = input$cd_p1, p2 = input$cd_p2,
          delta = input$cd_delta_binom
        ))
      },
      poisson = c(s, list(tipo = "poisson", l1 = input$cd_l1, l2 = input$cd_l2, delta = input$cd_delta_pois))
    )
  })

  servidor_contraste(input, output, session, "cd", cd_par, modelo_dos, c("n1", "n2"), m)

  output$cd_pob_plot <- renderPlot({
    m <- m()
    if (input$cd_dist == "normal") {
      s <- cd_par()
      mu <- c(s$mu1, s$mu2)
      sd <- c(s$sd1, s$sd2)
      nombres <- sprintf("Población %d: N(%s, %s)", 1:2, mu, sd)
      x <- seq(min(mu - 4 * sd), max(mu + 4 * sd), length.out = 500)
      d <- data.frame(
        x = c(x, x), y = c(dnorm(x, mu[1], sd[1]), dnorm(x, mu[2], sd[2])),
        pob = factor(rep(nombres, each = length(x)), levels = nombres)
      )
      g <- ggplot(d, aes(x, y, colour = pob, fill = pob)) +
        geom_area(position = "identity", alpha = 0.25) +
        geom_line(linewidth = 1) +
        labs(y = "Densidad")
    } else {
      if (input$cd_dist == "binomial") {
        ensayos <- input$cd_m
        p <- c(input$cd_p1, input$cd_p2)
        x <- 0:ensayos
        probs <- c(dbinom(x, ensayos, p[1]), dbinom(x, ensayos, p[2]))
        nombres <- sprintf("Población %d: B(%d, %s)", 1:2, ensayos, p)
      } else {
        l <- c(input$cd_l1, input$cd_l2)
        x <- 0:max(qpois(0.9995, l))
        probs <- c(dpois(x, l[1]), dpois(x, l[2]))
        nombres <- sprintf("Población %d: P(%s)", 1:2, l)
      }
      d <- data.frame(
        x = c(x, x), y = probs,
        pob = factor(rep(nombres, each = length(x)), levels = nombres)
      )
      g <- ggplot(d, aes(x, y, fill = pob)) +
        geom_col(position = position_dodge(width = 0.85), width = 0.8) +
        labs(y = "P(X = x)")
    }
    g +
      scale_fill_manual(values = setNames(m$series[c(7, 4)], nombres), name = NULL, aesthetics = c("colour", "fill")) +
      guides(fill = guide_legend(ncol = 1), colour = guide_legend(ncol = 1)) +
      labs(x = "x") +
      tema(m)
  })
}

shinyApp(ui, server)
