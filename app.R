# Simulaciones de procesos estadísticos aleatorios
# - Ley de los grandes números
# - Aproximación de la distribución binomial mediante la Poisson
# - Teorema central del límite
# - Intervalos de confianza para la media y para proporciones

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

# Recuadros de valores más compactos en pantallas pequeñas (< 576 px). En
# móvil bslib fuerza todas las rejillas a una columna; en las filas de valores
# se restablecen las 12 columnas para que los recuadros vayan de dos en dos.
# Los selectores llevan la clase .fila-valores para ganar en especificidad al
# CSS de bslib, que se carga después.
estilo_movil <- "
@media (max-width: 575.98px) {
  .bslib-grid.fila-valores {
    grid-template-columns: repeat(12, minmax(0, 1fr)) !important;
    gap: 0.5rem;
  }
  .fila-valores .bslib-value-box { margin-bottom: 0; }
  .fila-valores .bslib-value-box > .card-body { padding: 0; }
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
  "Aproximación binomial → Poisson",
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
}

shinyApp(ui, server)
