library(tidyverse)
library(lfe)
library(viridis)
library(gridExtra)
library(ggridges)
library(shiny)
library(shinydashboard)

ui <- bootstrapPage(dashboardBody(
  titlePanel("Econometrics Sandbox"),
  "Event Study Designs & Randomization Inference",
  helpText(paste0('Simulate randomization inference in event studies, with and without a pure control! ',
                  'Click "Simulate!" to generate data, and scroll down to see all the visualizations, for the data generating process, ',
                  'and for sliders to adjust fraction treated, true effects, and selection into treatment. ',
                  'Running randomization inference and generating the plots takes a few seconds per simulation.')),
  splitLayout(actionButton("runsim", "Simulate!")),
  fluidRow(box(plotOutput("esplot", height = "1000", width = "1000"))),
  withMathJax(),
  # SLIDERS!!!
  splitLayout(sliderInput("fraccntrl", "$$\\text{Fraction never treated } P$$", 0, 0.95, 0.5, step = 0.01),
              sliderInput("beta", "$$\\text{True effect } \\beta$$", 0, 10, 1, step = 0.1),
              sliderInput("beta_trend", "$$\\text{True effect growth } \\beta_{\\text{TREND}}$$", 0, 10, 0, step = 0.1)),
  splitLayout(sliderInput("level_w", "$$\\text{Selection on level } \\omega_{\\text{LEVEL}}$$", -1, 1, 0, step = 0.01),
              sliderInput("trend_w", "$$\\text{Selection on trend } \\omega_{\\text{TREND}}$$", -1, 1, 0, step = 0.01)),
  helpText('$$ \\begin{array}{lr} \\mu_{i} \\sim N(0, \\mathbf{I}_{4}) & d_{i} =
    \\mathbf{1} \\{ (1 - \\omega_{\\text{LEVEL}} - \\omega_{\\text{TREND}}) \\mu_{i0} +
    \\omega_{\\text{LEVEL}} \\mu_{i1} + \\omega_{\\text{TREND}} \\mu_{i2} > K_{1} \\} \\\\
    T_{i0} = \\mathbf{1} \\{ d_{i} = 1 \\; \\& \\; \\mu_{i3} < K_{2} \\} & T_{i2} = \\mathbf{1} \\{ d_{i} = 1 \\; \\& \\; \\mu_{i3} > K_{2} \\} \\\\
    T_{i} = 0 * T_{i0} + 2 * T_{i2} & d_{it} = d_{i} \\mathbf{1} \\{ t \\geq T_{i} \\} \\\\
    \\epsilon_{it,\\text{SHOCK}} \\sim N(0,1) & \\nu_{it,\\text{SHOCK}} = 0.5 \\nu_{it-1,\\text{SHOCK}} +
         \\epsilon_{it,\\text{SHOCK}} \\\\
    \\epsilon_{it,\\text{TREND}} \\sim N(0,1) & \\nu_{it,\\text{TREND}} = 0.5 \\nu_{it-1,\\text{TREND}} +
         \\epsilon_{it,\\text{TREND}}
    \\end{array} $$
  $$ y_{it} = (\\beta + (\\beta_{\\text{TREND}} (t - T_{i}) / T)) d_{it} + 0.4 \\mu_{i1} + 
         0.4 \\nu_{it,\\text{SHOCK}} +
         \\sum_{s = -T_{\\text{PRE}}}^{t} (0.4 \\mu_{i2} +
             0.4 \\nu_{is,\\text{TREND}}) / T $$'),
  helpText(paste0('\\(K_{1}\\) and \\(K_{2}\\) selected such that a share \\(P\\) of individuals are never treated, a share ',
                  '\\((1 - P)/2\\) are treated in period 0, and a share \\((1 - P)/2\\) are treated in period 2.'))
))

server <- function(input, output) {
  watercols5 <- viridis.map %>% filter(opt == "C") %>% `[`(,1:3) %>% sapply(function(x) { x[c(1, 64, 128, 192, 256)] }) %>%
    apply(MARGIN = 1, function(x) { paste0("#", round(x*255) %>% as.hexmode %>% paste(collapse = "")) })
  genplotsim <- eventReactive(input$runsim, {
    # input$fraccntrl
    # input$level_w
    # input$trend_w
    # input$beta
    # input$beta_trend
    ni <- 200
    fraccntrl <- input$fraccntrl %>% as.numeric
    if(!(fraccntrl >= 0 & fraccntrl <= 0.95)) stop("Please enter fraction control between 0 and 0.95.")
    ncntrl <- (round(ni*fraccntrl/2)*2) %>% as.integer
    ntreat0 <- ((ni - ncntrl)/2) %>% as.integer
    ntreat2 <- ntreat0
    if(ncntrl + ntreat0 + ntreat2 != ni) stop("Treated and control should sum to number of observations")
    
    t <- rep(-4:4, ni)
    nt <- 9
    i <- rep(1:ni, each = nt)
    mu0 <- rnorm(ni)
    mu1 <- rnorm(ni)
    mu2 <- rnorm(ni)
    level_w <- input$level_w %>% as.numeric
    trend_w <- input$trend_w %>% as.numeric
    # ASSIGNS BOTTOM ncntrl TO CONTROl, OTHERS TO TREATMENT
    di <- as.integer(rank((1 - level_w - trend_w) * mu0 + level_w * mu1 + trend_w * mu2) > ncntrl)
    # ASSIGNS FIRST HALF OF TREATED TO EARLY, LAST HALF TO LATE
    dilate <- as.integer(cumsum(di) <= ntreat0 & di)
    # MULTIPLIES BY NUMBER OF TIME PERIODS FOR EACH INDIVIDUAL
    mu1 <- mu1 %>% rep(each = nt)
    mu2 <- mu2 %>% rep(each = nt)
    di <- di %>% rep(each = nt)
    dilate <- dilate %>% rep(each = nt)
    
    # beta NUMBERS
    beta <- input$beta %>% as.numeric
    beta_trend <- input$beta_trend %>% as.numeric
    eps_shock <- rnorm(ni*nt)
    eps_trend <- rnorm(ni*nt)
    
    rho_shock <- 0.5
    rho_trend <- 0.5
    rhom_shock <- sapply(1:nt, function(s) {
      ifelse(s <= 1:nt, rho_shock^(1:nt - s), 0) * ifelse(s == 1, 1 / sqrt(1 - rho_shock^2), 1)
    }) * sqrt(1 - rho_shock^2)
    rhom_trend <- sapply(1:nt, function(s) {
      ifelse(s <= 1:nt, rho_trend^(1:nt - s), 0) * ifelse(s == 1, 1 / sqrt(1 - rho_trend^2), 1)
    }) * sqrt(1 - rho_trend^2)
    
    # _w NUMBERS
    y_level_w <- 0.4
    y_shock_w <- 0.4
    y_trend_w <- 0.4
    y_shocktrend_w  <- 0.4
    
    df <- data.frame(i = i, t = t, treat = di, treat2 = dilate,
                     mu1 = mu1, mu2 = mu2, eps_shock = eps_shock, eps_trend = eps_trend) %>%
      group_by(i) %>% arrange(t) %>%
      mutate(nu_shock = rhom_shock %*% eps_shock,
             nu_trend = rhom_trend %*% eps_trend,
             tau = ifelse(rep(all(treat == 0), n()),
                          rep(-1, n()),
                          t - 2*treat2)) %>%
      mutate(y = (beta + (beta_trend / nt)*tau) * (tau >= 0) + y_shock_w * nu_shock + y_level_w * mu1 +
               cumsum(y_trend_w * mu2 + y_shocktrend_w * nu_trend) / nt) %>%
      ungroup %>%
      # TRIM tau
      mutate(tau = case_when(tau > 2 ~ 2, tau < -4 ~ -4, T ~ tau)) %>%
      # FACTOR
      mutate(tau = factor(tau, levels = c(-1, unique(tau) %>% subset(. != -1)))) %>%
      arrange(treat, treat2, i, t)
    
    rawdata <- df %>% group_by(treat, treat2, t) %>%
      summarise(se = sd(y) / sqrt(n()), y = mean(y)) %>%
      mutate(group = case_when(treat == 0 ~ "Control", treat2 == 0 ~ "Treated 0", T ~ "Treated 2"))
    regeventstudy <- felm(y ~ tau | i + t | 0 | i, data = df) %>% summary %>% `$`("coef")
    regeventstudydf <- data.frame(tau = row.names(regeventstudy) %>% gsub(pattern = "tau", replacement = "") %>% as.numeric,
                                  coef = regeventstudy[,"Estimate"], se = regeventstudy[,"Cluster s.e."]) %>%
      mutate(lci = coef + se * qnorm(.025), uci = coef + se * qnorm(.975)) %>%
      bind_rows(data.frame(tau = -1, coef = 0))
    regeventstudync <- felm(y ~ tau | i + t | 0 | i, data = df %>% filter(treat == 1)) %>% summary %>% `$`("coef")
    regeventstudyncdf <- data.frame(tau = row.names(regeventstudync) %>% gsub(pattern = "tau", replacement = "") %>% as.numeric,
                                    coef = regeventstudync[,"Estimate"], se = regeventstudync[,"Cluster s.e."]) %>%
      mutate(lci = coef + se * qnorm(.025), uci = coef + se * qnorm(.975)) %>%
      bind_rows(data.frame(tau = -1, coef = 0))
    
    # RANDOMIZATION INFERENCE MATRIX GENERATED ONCE !
    rimatrix <- rnorm(1e3 * ni) %>% matrix(nrow = 1e3, ncol = ni)
    for(sim in 1:1e3) {
      rimatrix[sim,] <- rank(rimatrix[sim,])
    }
    
    # DO DEMEANING BY FIXED EFFECTS OF tau BEFORE PERMUTING --> DON'T NEED TO APPLY felm!
    dfdm <- df
    rhs <- paste0("`tau", levels(df$tau) %>% `[`(2:length(.)), "`")
    for(k in 2:length(levels(df$tau))) {
      tauk <- levels(dfdm$tau)[k]
      dfdm[,paste0("tau", tauk)] <- dfdm$tau == tauk
      dfdm[,paste0("tau", tauk)] <- felm(as.formula(paste0("`tau", tauk, "` ~ 1 | i + t")), data = dfdm)$resid %>% as.numeric
    }
    # GENERATE X MATRIX (SO CAN DO RANDOMIZATION INFERENCE BY REORDERING y'S)
    xmatrix <- dfdm[,rhs %>% gsub(pattern = "`", replacement = "")] %>% as.matrix
    xpxinvxp <- solve(t(xmatrix) %*% xmatrix) %*% t(xmatrix)
    
    dfdmnc <- df %>% subset(treat == 1)
    rhs <- paste0("`tau", levels(df$tau) %>% `[`(2:length(.)), "`")
    for(k in 2:length(levels(df$tau))) {
      tauk <- levels(dfdmnc$tau)[k]
      dfdmnc[,paste0("tau", tauk)] <- dfdmnc$tau == tauk
      dfdmnc[,paste0("tau", tauk)] <- felm(as.formula(paste0("`tau", tauk, "` ~ 1 | i + t")), data = dfdmnc)$resid %>% as.numeric
    }
    xmatrixnc <- dfdmnc[,rhs %>% gsub(pattern = "`", replacement = "")] %>% as.matrix
    xpxinvxpnc <- solve(t(xmatrixnc) %*% xmatrixnc) %*% t(xmatrixnc)
    
    ri1 <- lapply(1:1000, function(sim) {
      # RANDOMIZED REASSIGNMENT OF tau
      isim <- (rimatrix[sim,] %>% rep(each = nt)) * nt + ((1:nt) - nt)
      coefsim <- xpxinvxp[,isim] %*% as.numeric(dfdm$y)
      regsimdf <- data.frame(tau = row.names(coefsim) %>% gsub(pattern = "tau", replacement = "") %>% 
                               gsub(pattern = "`", replacement = "") %>% as.numeric,
                             coef = coefsim %>% as.numeric %>% unname, sim = sim)
    }) %>% bind_rows()
    ri1df <- ri1 %>% group_by(tau) %>%
      summarise(lci = quantile(coef, .975), uci = quantile(coef, .025)) %>%
      ungroup %>%
      left_join(regeventstudydf %>% select(tau, coef)) %>%
      mutate(lci = coef - lci, uci = coef - uci)
    
    ri2 <- lapply(1:1000, function(sim) {
      # RANDOMIZES PRESERVING ASSIGNMENT TO CONTROL
      isim <- c(rimatrix[sim,] %>% subset(. > ncntrl))
      if(ncntrl > 0) isim <- c(1:ncntrl, isim)
      isim <- (isim %>% rep(each = nt)) * nt + ((1:nt) - nt)
      coefsim <- xpxinvxp[,isim] %*% as.numeric(dfdm$y)
      regsimdf <- data.frame(tau = row.names(coefsim) %>% gsub(pattern = "tau", replacement = "") %>% 
                               gsub(pattern = "`", replacement = "") %>% as.numeric,
                             coef = coefsim %>% as.numeric %>% unname, sim = sim)
    }) %>% bind_rows()
    ri2df <- ri2 %>% group_by(tau) %>%
      summarise(lci = quantile(coef, .975), uci = quantile(coef, .025)) %>%
      ungroup %>%
      left_join(regeventstudydf %>% select(tau, coef)) %>%
      mutate(lci = coef - lci, uci = coef - uci)
    
    rinc <- lapply(1:1000, function(sim) {
      isim <- rimatrix[sim,] %>% subset(. > ncntrl) - ncntrl
      isim <- (isim %>% rep(each = nt)) * nt + ((1:nt) - nt)
      coefsim <- xpxinvxpnc[,isim] %*% as.numeric(dfdmnc$y)
      regsimdf <- data.frame(tau = row.names(coefsim) %>% gsub(pattern = "tau", replacement = "") %>% 
                               gsub(pattern = "`", replacement = "") %>% as.numeric,
                             coef = coefsim %>% as.numeric %>% unname, sim = sim)
    }) %>% bind_rows()
    rincdf <- rinc %>% group_by(tau) %>%
      summarise(lci = quantile(coef, .975), uci = quantile(coef, .025)) %>%
      ungroup %>%
      left_join(regeventstudyncdf %>% select(tau, coef)) %>%
      mutate(lci = coef - lci, uci = coef - uci)
    
    coefrange <- range(c(ri1$coef, ri2$coef, rinc$coef, regeventstudydf$lci, regeventstudydf$uci,
                         regeventstudyncdf$lci, regeventstudyncdf$uci), na.rm = T)
    plotbase_size = 16
    g1ridge <- ggplot() +
      geom_density_ridges(data = bind_rows(ri1 %>% mutate(group = paste0("1", tau), model = "1treated"),
                                           ri2 %>% mutate(group = paste0("2", tau), model = "2timetreated")),
                          aes(x = coef, y = tau, group = group, fill = model), alpha = 0.5) +
      geom_point(data = regeventstudydf, aes(x = coef, y = tau), size = 6) +
      scale_y_continuous(breaks = -4:2, labels = -4:2, name = "Event time") +
      xlab("Point estimate") + xlim(coefrange[1], coefrange[2]) +
      scale_fill_manual(labels = c("Randomize\nTreated and\nTime Treated", "Randomize\nTime Treated"),
                        values = c("1treated" = watercols5[2], "2timetreated" = watercols5[3]),
                        name = NULL) +
      theme_bw(base_size = plotbase_size) + theme(legend.position = "top", legend.direction = "horizontal") +
      ggtitle("Event Study Randomization Inference", subtitle = "Include never treated")
    g2ridge <- ggplot() +
      geom_density_ridges(data = rinc %>% mutate(group = paste0("3", tau), model = "3timetreated"),
                          aes(x = coef, y = tau, group = group, fill = model), alpha = 0.5) +
      geom_point(data = regeventstudyncdf, aes(x = coef, y = tau), size = 6) +
      scale_y_continuous(breaks = -4:2, labels = -4:2, name = "Event time") +
      xlab("Point estimate") + xlim(coefrange[1], coefrange[2]) +
      scale_fill_manual(labels = c("Randomize\nTime Treated"),
                        values = c("3timetreated" = watercols5[3]),
                        name = NULL) +
      theme_bw(base_size = plotbase_size) + theme(legend.position = "top", legend.direction = "horizontal") +
      ggtitle("Event Study Randomization Inference", subtitle = "Drop never treated")
    grawdata <- ggplot() +
      geom_point(data = rawdata, aes(x = t, y = y, col = group), size = 4) +
      geom_line(data = rawdata, aes(x = t, y = y, col = group), size = 2) +
      geom_vline(xintercept = -0.5, col = watercols5[4], linetype = 2, size = 2) +
      geom_vline(xintercept = 1.5, col = watercols5[5], linetype = 2, size = 2) +
      scale_color_manual(values = c("Control" = watercols5[1], "Treated 0" = watercols5[4], "Treated 2" = watercols5[5]),
                         name = NULL) +
      xlab("Time") + ylab(NULL) +
      theme_bw(base_size = plotbase_size)+ theme(legend.position = "top", legend.direction = "horizontal") + ggtitle("Raw Means")
    g1event <- ggplot() +
      geom_point(data = regeventstudydf, aes(x = tau, y = coef), size = 6) +
      geom_errorbar(data = regeventstudydf, aes(x = tau, ymin = lci, ymax = uci), size = 3, width = 0.3) +
      scale_x_continuous(breaks = -4:2, labels = -4:2, name = "Event time") +
      ylab("Point estimate") + ylim(coefrange[1], coefrange[2]) +
      theme_bw(base_size = plotbase_size) + ggtitle("Event Study", subtitle = "Include never treated")
    g2event <- ggplot() +
      geom_point(data = regeventstudyncdf, aes(x = tau, y = coef), size = 6) +
      geom_errorbar(data = regeventstudyncdf, aes(x = tau, ymin = lci, ymax = uci), size = 3, width = 0.3) +
      scale_x_continuous(breaks = -4:2, labels = -4:2, name = "Event time") +
      ylab("Point estimate") + ylim(coefrange[1], coefrange[2]) +
      theme_bw(base_size = plotbase_size) + ggtitle("Event Study", subtitle = "Drop never treated")
    grid.arrange(
      grobs = list(grawdata, g1event, g1ridge, g2event, g2ridge),
      layout_matrix = rbind(c(1, 1), c(2, 3), c(4, 5))
    )
  })
  
  output$esplot <- renderPlot({
    genplotsim()
  })
}

shinyApp(ui, server)

