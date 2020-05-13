library(tidyverse)
library(lfe)
library(viridis)
library(shiny)
library(shinydashboard)

ui <- bootstrapPage(dashboardBody(
  titlePanel("Econometrics Sandbox"),
  "Event Study Designs & Co.",
  fluidRow(box(plotOutput("esplot", height = "500", width = "1000"))),
  helpText(paste0('Simulate event studies, difference-in-differences, and regression discontinuity in time! ',
                  'Click "Simulate!" to generate data, and scroll down for the data generating process ',
                  'and additional sliders to adjust sample size, true effects, serial correlation, selection, and noise.')),
  splitLayout(actionButton("runsim", "Simulate!")),
  withMathJax(),
  # SLIDERS!!!
  splitLayout(sliderInput("nt", "$$\\text{# of time periods } T$$", 2, 100, 20, step = 1),
              # textInput("nt", "$$\\text{# of time periods } T$$", 20),
              sliderInput("ntpre", "$$\\text{# of pre-periods } T_{\\text{PRE}}$$", 1, 99, 10, step = 1),
              # textInput("ntpre", "$$\\text{# of pre-periods } T_{\\text{PRE}}$$", 10),
              sliderInput("beta", "$$\\text{True effect } \\beta$$", 0, 10, 1, step = 0.1),
              # textInput("beta", "$$\\text{True effect } \\beta$$", 1),
              sliderInput("beta_trend", "$$\\text{True effect growth } \\beta_{\\text{TREND}}$$", 0, 10, 2, step = 0.1)),
              # textInput("beta_trend", "$$\\text{True effect growth } \\beta_{\\text{TREND}}$$", 2)),
  splitLayout(sliderInput("ni", "$$\\text{# of individuals } N$$", 50, 5000, 1000, step = 50),
              # textInput("ni", "$$\\text{# of individuals } N$$", 1000),
              sliderInput("rho_shock", "$$\\text{Serial correlation error } \\rho_{\\text{SHOCK}}$$", -0.95, 0.95, 0.9),
              # textInput("rho_shock", "$$\\text{Serial correlation error } \\rho_{\\text{SHOCK}}$$", 0.9),
              sliderInput("rho_trend", "$$\\text{Serial correlation trend error } \\rho_{\\text{TREND}}$$", -0.95, 0.95, 0.9),
              # textInput("rho_trend", "$$\\text{Serial correlation trend error } \\rho_{\\text{TREND}}$$", 0.9),
              sliderInput("hrdit", "$$\\text{RDiT bandwidth}$$", 1, 50, 3, step = 1)),
              # textInput("hrdit", "$$\\text{RDiT bandwidth}$$", 3)),
  splitLayout(sliderInput("level_w", "$$\\text{Selection on level } \\omega_{\\text{LEVEL}}$$", -1, 1, 0, step = 0.01),
              # textInput("level_w", "$$\\text{Selection on level } \\omega_{\\text{LEVEL}}$$", 0),
              sliderInput("trend_w", "$$\\text{Selection on trend } \\omega_{\\text{TREND}}$$", -1, 1, 0, step = 0.01),
              # textInput("trend_w", "$$\\text{Selection on trend } \\omega_{\\text{TREND}}$$", 0),
              sliderInput("shock_w", "$$\\text{Selection on shock } \\omega_{\\text{SHOCK}}$$", -1, 1, 0, step = 0.01),
              # textInput("shock_w", "$$\\text{Selection on shock } \\omega_{\\text{SHOCK}}$$", 0),
              sliderInput("shocktrend_w", "$$\\text{Selection on trend shock } \\omega_{\\text{TRENDSHOCK}}$$", -1, 1, 0, step = 0.01)),
              # textInput("shocktrend_w", "$$\\text{Selection on trend shock } \\omega_{\\text{TRENDSHOCK}}$$", 0)),
  splitLayout(sliderInput("y_level_w", "$$\\text{Variation in level } \\lambda_{\\text{LEVEL}}$$", 0, 10, 1, step = 0.1),
              # textInput("y_level_w", "$$\\text{Variation in level } \\lambda_{\\text{LEVEL}}$$", 1),
              sliderInput("y_trend_w", "$$\\text{Variation in trend } \\lambda_{\\text{TREND}}$$", 0, 10, 1, step = 0.1),
              # textInput("y_trend_w", "$$\\text{Variation in trend } \\lambda_{\\text{TREND}}$$", 1),
              sliderInput("y_shock_w", "$$\\text{Variation in shock } \\lambda_{\\text{SHOCK}}$$", 0, 10, 1, step = 0.1),
              # textInput("y_shock_w", "$$\\text{Variation in shock } \\lambda_{\\text{SHOCK}}$$", 1),
              sliderInput("y_shocktrend_w", "$$\\text{Variation in trend shock } \\lambda_{\\text{TRENDSHOCK}}$$", 0, 10, 1, step = 0.1)),
              # textInput("y_shocktrend_w", "$$\\text{Variation in trend shock } \\lambda_{\\text{TRENDSHOCK}}$$", 1)),
  helpText('$$ \\begin{array}{lr} \\mu_{i} \\sim N(0, \\mathbf{I}_{5}) & d_{it} =
    \\mathbf{1} \\{ t \\geq 0 \\; \\& \\; (1 - \\omega_{\\text{LEVEL}} - \\omega_{\\text{TREND}} -
    \\omega_{\\text{SHOCK}} - \\omega_{\\text{TRENDSHOCK}}) \\mu_{i0} \\\\
    & + \\omega_{\\text{LEVEL}} \\mu_{i1} +
    \\omega_{\\text{TREND}} \\mu_{i2} + \\omega_{\\text{SHOCK}} \\mu_{i3} + \\omega_{\\text{TRENDSHOCK}} \\mu_{i4} > 0 \\} \\\\
    \\epsilon_{it,\\text{SHOCK}} \\sim N(0,1) & \\nu_{it,\\text{SHOCK}} = \\rho_{\\text{SHOCK}} \\nu_{it-1,\\text{SHOCK}} +
         \\mathbf{1} \\{ t \\neq 0 \\} \\epsilon_{it,\\text{SHOCK}} +
         \\mathbf{1} \\{ t = 0 \\} \\mu_{i3} \\\\
    \\epsilon_{it,\\text{TREND}} \\sim N(0,1) & \\nu_{it,\\text{TREND}} = \\rho_{\\text{TREND}} \\nu_{it-1,\\text{TREND}} +
         \\mathbf{1} \\{ t \\neq 0 \\} \\epsilon_{it,\\text{TREND}} +
         \\mathbf{1} \\{ t = 0 \\} \\mu_{i4}
    \\end{array} $$
  $$ y_{it} = (\\beta + (\\beta_{\\text{TREND}} t / T)) d_{it} + \\lambda_{\\text{LEVEL}} \\mu_{i1} + 
         \\lambda_{\\text{SHOCK}} \\nu_{it,\\text{SHOCK}} +
         \\sum_{s = -T_{\\text{PRE}}}^{t} (\\lambda_{\\text{TREND}} \\mu_{i2} +
             \\lambda_{\\text{SHOCKTREND}} \\nu_{is,\\text{TREND}}) / T $$')
))

server <- function(input, output) {
  genplotsim <- eventReactive(input$runsim, {
    watercols5 <- viridis.map %>% filter(opt == "C") %>% `[`(,1:3) %>% sapply(function(x) { x[c(1, 64, 128, 192, 256)] }) %>%
      apply(MARGIN = 1, function(x) { paste0("#", round(x*255) %>% as.hexmode %>% paste(collapse = "")) })
    
    # ni = NUMBER OF INDIVIDUALS, nt = NUMBER OF TIME PERIODS
    ni <- input$ni %>% as.numeric
    if(!(ni >= 20 & ni <= 1e4)) stop("Please enter at least 20 and at most 10000 individuals.")
    nt <- input$nt %>% as.numeric
    if(!(nt >= 2 & nt <= 1e3)) stop("Please enter at least 2 and at most 1000 time periods")
    ntpre <- input$ntpre %>% as.numeric
    if(!(ntpre >= 1 & ntpre <= nt - 1)) stop("Please enter at least 1 and at most T - 1 pre-periods")
    # t INDEXES TIME, i INDEXES INDIVIDUALS
    t <- rep(1:nt, ni) - (ntpre + 1)
    i <- rep(1:ni, each = nt)
    # RD BANDWIDTH
    hrdit <- input$hrdit %>% as.numeric
    if(!(hrdit >= 1 & hrdit <= min(max(abs(t)) + 1, max(abs(-t))))) {
      stop("Please enter RDiT bandwidth at least 1 and at most the min of the number of pre-periods and the number of post-periods")
    }
    eta <- -0.5 + (t >= 0)
    mu0 <- rep(rnorm(ni), each = nt)
    mu1 <- rep(rnorm(ni), each = nt)
    mu2 <- rep(rnorm(ni), each = nt)
    mu3 <- rep(rnorm(ni), each = nt)
    mu4 <- rep(rnorm(ni), each = nt)
    # SUM OF _w <= 1
    level_w <- input$level_w %>% as.numeric
    trend_w <- input$trend_w %>% as.numeric
    shock_w <- input$shock_w %>% as.numeric
    shocktrend_w <- input$shocktrend_w %>% as.numeric
    # if(!(level_w >= 0 & trend_w >= 0 & shock_w >= 0 & shocktrend_w >= 0 & (level_w + trend_w + shock_w + shocktrend_w) <= 1)) {
    #   stop("Selection coefficients must...")
    # }
    d <- as.integer(eta > 0 & ((1 - level_w - trend_w - shock_w - shocktrend_w) * mu0 +
                                 level_w * mu1 + trend_w * mu2 + shock_w * mu3 + shocktrend_w * mu4) > 0)
    # beta NUMBERS
    beta <- input$beta %>% as.numeric
    beta_trend <- input$beta_trend %>% as.numeric
    eps_shock <- ifelse(t == 0, mu3, rnorm(ni*nt))
    eps_trend <- ifelse(t == 0, mu4, rnorm(ni*nt))
    # -1 TO 1 (CAN'T BE -1 OR 1)
    rho_shock <- input$rho_shock %>% as.numeric
    if(!(rho_shock > -1 & rho_shock < 1)) stop("Please enter serial correlation in errors strictly between -1 and 1")
    # -1 TO 1 (CAN'T BE -1 OR 1)
    rho_trend <- input$rho_trend %>% as.numeric
    if(!(rho_trend > -1 & rho_trend < 1)) stop("Please enter serial correlation in trend errors strictly between -1 and 1")
    rhom_shock <- sapply(1:nt, function(i) {
      ifelse(i <= 1:nt, rho_shock^(1:nt - i), 0) * ifelse(i == 1, 1 / sqrt(1 - rho_shock^2), 1)
    }) * sqrt(1 - rho_shock^2)
    rhom_trend <- sapply(1:nt, function(i) {
      ifelse(i <= 1:nt, rho_trend^(1:nt - i), 0) * ifelse(i == 1, 1 / sqrt(1 - rho_trend^2), 1)
    }) * sqrt(1 - rho_trend^2)
    # _w NUMBERS
    y_level_w <- input$y_level_w %>% as.numeric
    y_shock_w <- input$y_shock_w %>% as.numeric
    y_trend_w <- input$y_trend_w %>% as.numeric
    y_shocktrend_w  <- input$y_shocktrend_w %>% as.numeric
    df <- data.frame(i = i, t = t, d = d, eps_shock = eps_shock, eps_trend = eps_trend,
                     mu0 = mu0, mu1 = mu1, mu2 = mu2, mu3 = mu3, mu4 = mu4) %>%
      group_by(i) %>% arrange(t) %>%
      mutate(nu_shock = rhom_shock %*% eps_shock,
             nu_trend = rhom_trend %*% eps_trend,
             tau = ifelse(rep(all(d == 0), n()),
                          rep(-1, n()),
                          t - min(ifelse(d == 0, Inf, t)))) %>%
      mutate(y = (beta + (beta_trend / nt)*tau) * d + y_shock_w * nu_shock + y_level_w * mu1 +
               cumsum(y_trend_w * mu2 + y_shocktrend_w * nu_trend) / nt) %>%
      ungroup %>%
      mutate(tau = factor(tau, levels = c(-1, unique(tau) %>% subset(. != -1))))
    
    regeventstudy <- felm(y ~ tau | i + t | 0 | i, data = df) %>% summary %>% `$`("coef")
    regeventstudydf <- data.frame(tau = row.names(regeventstudy) %>% gsub(pattern = "tau", replacement = "") %>% as.numeric,
                                  coef = regeventstudy[,"Estimate"], se = regeventstudy[,"Cluster s.e."]) %>%
      mutate(lci = coef + se * qnorm(.025), uci = coef + se * qnorm(.975)) %>%
      bind_rows(data.frame(tau = -1, coef = 0))
    regdd <- felm(y ~ d | i + t | 0 | i, data = df) %>% summary %>% `$`("coef")
    regdddf <- data.frame(coef = regdd["d","Estimate"], se = regdd["d", "Cluster s.e."]) %>%
      mutate(lci = coef + se * qnorm(.025), uci = coef + se * qnorm(.975))
    regrdit <- felm(y ~ d | i + t | 0 | i, data = df %>%
                      filter(as.numeric(as.character(tau)) >= -hrdit & as.numeric(as.character(tau)) < hrdit)) %>%
      summary %>% `$`("coef")
    regrditdf <- data.frame(coef = regrdit["d","Estimate"], se = regrdit["d", "Cluster s.e."]) %>%
      mutate(lci = coef + se * qnorm(.025), uci = coef + se * qnorm(.975))
    
    ggplot() +
      # DID
      annotate(geom = "segment", x = min(regeventstudydf$tau), xend = 0, y = 0, yend = 0, col = watercols5[2], size = 1) +
      geom_segment(data = regdddf, aes(y = coef, yend = coef, col = "DID"), x = 0, xend = max(regeventstudydf$tau), size = 1) +
      annotate(geom = "segment", x = 0, xend = max(regeventstudydf$tau), y = regdddf$lci, yend = regdddf$lci, col = watercols5[2],
               linetype = 2, size = 1) +
      annotate(geom = "segment", x = 0, xend = max(regeventstudydf$tau), y = regdddf$uci, yend = regdddf$uci, col = watercols5[2],
               linetype = 2, size = 1) +
      # RDIT
      annotate(geom = "segment", x = pmax(min(regeventstudydf$tau), -hrdit), xend = 0, y = 0, yend = 0, col = watercols5[3], size = 1) +
      geom_segment(data = regrditdf, aes(y = coef, yend = coef, col = "RDiT"), x = 0, xend = pmin(max(regeventstudydf$tau), hrdit - 1),
                   size = 1) +
      annotate(geom = "segment", x = 0, xend = pmin(max(regeventstudydf$tau), hrdit - 1),
               y = regrditdf$lci, yend = regrditdf$lci, col = watercols5[3], linetype = 2, size = 1) +
      annotate(geom = "segment", x = 0, xend = pmin(max(regeventstudydf$tau), hrdit - 1),
               y = regrditdf$uci, yend = regrditdf$uci, col = watercols5[3], linetype = 2, size = 1) +
      # EVENT STUDY
      geom_point(data = regeventstudydf, aes(x = tau, y = coef, col = "Event study"), size = 3) +
      geom_linerange(data = regeventstudydf, aes(x = tau, ymin = lci, ymax = uci), size = 1) +
      theme_bw() +
      scale_colour_manual(name = NULL, values=c("Event study" = "black", "DID" = watercols5[2], "RDiT" = watercols5[3])) +
      theme(legend.position = "top", legend.direction = "horizontal",
            legend.text = element_text(size = 16), axis.text = element_text(size = 16),
            axis.title = element_text(size = 16)) +
      xlab("Event time") + ylab(NULL)
  })
  
  output$esplot <- renderPlot({
    genplotsim()
  })
}

shinyApp(ui, server)

