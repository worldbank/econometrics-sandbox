library(tidyverse)
library(lfe)
library(viridis)
library(shiny)
library(shinydashboard)

# input <- list(nobs = 1e2,
#               nsims = 4e2,
#               sigeps00 = 0.2,
#               sigeps10 = 0.2,
#               sigeps11 = 0.2,
#               # rho0010 = 0,
#               # rho0011 = 0,
#               # rho1011 = 0,
#               beta1 = 0,
#               beta2 = 0,
#               beta3 = 0,
#               beta4 = 0,
#               beta5 = 0,
#               beta6 = -1,
#               tau = 1,
#               pattrit0 = 0.1,
#               pattrit1 = 0.4,
#               rhoattrit = 0)


ui <- bootstrapPage(dashboardBody(
  titlePanel("Econometrics Sandbox"),
  "Attrition",
  fluidRow(box(plotOutput("esplot", height = "500", width = "1000"))),
  helpText(paste0('Simulate attrition in an experiment when running difference-in-differences! ',
                  'Click "Simulate!" to generate data, and scroll down for the data generating process ',
                  'and additional sliders to adjust sample size, noise, attrition, and selection into attrition.')),
  splitLayout(actionButton("runsim", "Simulate!")),
  withMathJax(),
  # SLIDERS!!!
  splitLayout(sliderInput("nobs", "$$\\text{# of individuals}$$", 20, 200, 100, step = 2),
              sliderInput("nsims", "$$\\text{# of simulations}$$", 100, 1000, 400, step = 10),
              sliderInput("tau", "$$\\text{True effect } \\tau$$", 0, 10, 1, step = 0.1),
              sliderInput("sigeps00", "$$\\text{SD of baseline outcome } \\sigma_{00}$$", 0.01, 1, 0.2, step = 0.01)),
  splitLayout(sliderInput("pattrit0", "$$\\text{Control attrition } P(0)$$", 0, 0.99, 0.2, step = 0.01),
              sliderInput("pattrit1", "$$\\text{Treatment attrition } P(1)$$", 0, 0.99, 0.2, step = 0.01),
              sliderInput("sigeps10", "$$\\text{SD of control outcome } \\sigma_{10}$$", 0.01, 1, 0.2, step = 0.01),
              sliderInput("sigeps11", "$$\\text{SD of treated outcome } \\sigma_{11}$$", 0.01, 1, 0.2, step = 0.01)),
  splitLayout(sliderInput("beta1", "$$\\text{Baseline outcome of control attriters } \\beta_{1}$$", -5, 5, 0, step = 0.1),
              sliderInput("beta2", "$$\\text{Baseline outcome of treatment attriters } \\beta_{2}$$", -5, 5, 0, step = 0.1),
              sliderInput("beta3", "$$\\text{Post outcome of control attriters } \\beta_{3}$$", -5, 5, 0, step = 0.1),
              sliderInput("beta4", "$$\\text{Post outcome of treatment attriters } \\beta_{4}$$", -5, 5, 0, step = 0.1)),
  splitLayout(sliderInput("rhoattrit", "$$\\text{Correlation in attrition selection } \\rho$$", -1, 1, 1, step = 0.01),
              "",
              "",
              ""),
  # splitLayout(sliderInput("beta3", "$$\\text{Control outcome of control attriters } \\beta_{3}$$", -5, 5, 0, step = 0.1),
  #             sliderInput("beta4", "$$\\text{Control outcome of treatment attriters } \\beta_{4}$$", -5, 5, 0, step = 0.1),
  #             sliderInput("beta5", "$$\\text{Treated outcome of control attriters } \\beta_{5}$$", -5, 5, 0, step = 0.1),
  #             sliderInput("beta6", "$$\\text{Treated outcome of treatment attriters } \\beta_{6}$$", -5, 5, 0, step = 0.1)),
  helpText(paste0('All individuals \\(i\\) observed at baseline with outcome \\(Y_{i0}\\). ',
                  'Individuals who do not attrit (\\(S_{i} = 0\\)) have post treatment outcome ',
                  '\\(Y_{i1} = D_{i} Y_{i1}(1) + (1 - D_{i}) Y_{i1}(0)\\) observed. Attrition ',
                  '\\(S_{i} = D_{i} S_{i}(1) + (1 - D_{i}) S_{i}(0)\\). "OLS" estimates difference-in-differences specification ',
                  '\\(Y_{it} = \\alpha + \\delta D_{i} + \\gamma t + \\tau D_{i} * t + \\epsilon_{it}\\). "Drop attriters" estimates the ',
                  'same specification but drops all observations that attrit in period 1; this is equivalent to including individual ',
                  'fixed effects. Lee (2009) bounds implement Lee bounds ',
                  'using \\(Y_{i1} - Y_{i0}\\) as the outcome. ',
                  '$$ \\begin{array}{lr} Y_{i0} = \\epsilon_{i0} + \\beta_{1} s_{i}(0) + \\beta_{2} s_{i}(1) &
                  Y_{i1}(0) = \\epsilon_{i1}(0) + \\beta_{3} s_{i}(0) + \\beta_{4} s_{i}(1) \\\\
                  Y_{i1}(1) = \\epsilon_{i1}(1) + \\beta_{3} s_{i}(0) + \\beta_{4} s_{i}(1) + \\tau & \\\\
                  S_{i}(0) = \\mathbf{1} \\{ s_{i}(0) > \\Phi^{-1}(1 - P(0)) \\} &
                  S_{i}(1) = \\mathbf{1} \\{ s_{i}(1) > \\Phi^{-1}(1 - P(1)) \\} \\\\
                  (\\epsilon_{i0}, \\epsilon_{i1}(0), \\epsilon_{i1}(1)) \\sim N \\left( \\mathbf{0},
                  \\begin{array}{ccc} \\sigma_{00}^{2} & 0 & 0 \\\\ 0 & \\sigma_{10}^{2} & 0 \\\\ 0 & 0 & \\sigma_{11}^{2} \\end{array}
                  \\right) & (s_{i}(0), s_{i}(1)) \\sim N \\left( \\mathbf{0},
                  \\begin{array}{cc} 1 & \\rho \\\\ \\rho & 1 \\end{array}
                  \\right) \\end{array} $$'))
  # helpText('$$ \\begin{array}{lr} \\mu_{i} \\sim N(0, \\mathbf{I}_{5}) & d_{it} =
  #   \\mathbf{1} \\{ t \\geq 0 \\; \\& \\; (1 - \\omega_{\\text{LEVEL}} - \\omega_{\\text{TREND}} -
  #   \\omega_{\\text{SHOCK}} - \\omega_{\\text{TRENDSHOCK}}) \\mu_{i0} \\\\
  #   & + \\omega_{\\text{LEVEL}} \\mu_{i1} +
  #   \\omega_{\\text{TREND}} \\mu_{i2} + \\omega_{\\text{SHOCK}} \\mu_{i3} + \\omega_{\\text{TRENDSHOCK}} \\mu_{i4} > 0 \\} \\\\
  #   \\epsilon_{it,\\text{SHOCK}} \\sim N(0,1) & \\nu_{it,\\text{SHOCK}} = \\rho_{\\text{SHOCK}} \\nu_{it-1,\\text{SHOCK}} +
  #        \\mathbf{1} \\{ t \\neq 0 \\} \\epsilon_{it,\\text{SHOCK}} +
  #        \\mathbf{1} \\{ t = 0 \\} \\mu_{i3} \\\\
  #   \\epsilon_{it,\\text{TREND}} \\sim N(0,1) & \\nu_{it,\\text{TREND}} = \\rho_{\\text{TREND}} \\nu_{it-1,\\text{TREND}} +
  #        \\mathbf{1} \\{ t \\neq 0 \\} \\epsilon_{it,\\text{TREND}} +
  #        \\mathbf{1} \\{ t = 0 \\} \\mu_{i4}
  #   \\end{array} $$
  # $$ y_{it} = (\\beta + (\\beta_{\\text{TREND}} t / T)) d_{it} + \\lambda_{\\text{LEVEL}} \\mu_{i1} +
  #        \\lambda_{\\text{SHOCK}} \\nu_{it,\\text{SHOCK}} +
  #        \\sum_{s = -T_{\\text{PRE}}}^{t} (\\lambda_{\\text{TREND}} \\mu_{i2} +
  #            \\lambda_{\\text{SHOCKTREND}} \\nu_{is,\\text{TREND}}) / T $$')
))

server <- function(input, output) {
  genplotsim <- eventReactive(input$runsim, {
    nobs <- input$nobs %>% as.numeric
    nsims <- input$nsims %>% as.numeric
    sigeps00 <- input$sigeps00 %>% as.numeric
    sigeps10 <- input$sigeps10 %>% as.numeric
    sigeps11 <- input$sigeps11 %>% as.numeric
    # rho0010 <- input$rho0010 %>% as.numeric
    # rho0011 <- input$rho0011 %>% as.numeric
    # rho1011 <- input$rho1011 %>% as.numeric
    beta1 <- input$beta1 %>% as.numeric
    beta2 <- input$beta2 %>% as.numeric
    beta3 <- input$beta3 %>% as.numeric
    beta4 <- input$beta4 %>% as.numeric
    tau <- input$tau %>% as.numeric
    pattrit0 <- input$pattrit0 %>% as.numeric
    pattrit1 <- input$pattrit1 %>% as.numeric
    rhoattrit <- input$rhoattrit %>% as.numeric
    
    errors <- (rnorm(nsims * nobs * 5) %>% matrix(ncol = 5)) %*%
      # chol(matrix(c(sigeps00^2, rho0010 * sigeps00 * sigeps10 - .Machine$double.eps, rho0011 * sigeps00 * sigeps11 - .Machine$double.eps, 0, 0,
      #               rho0010 * sigeps00 * sigeps10 - .Machine$double.eps, sigeps10^2, rho1011 * sigeps10 * sigeps11 - .Machine$double.eps, 0, 0,
      #               rho0011 * sigeps00 * sigeps11 - .Machine$double.eps, rho1011 * sigeps10 * sigeps11 - .Machine$double.eps, sigeps11^2, 0, 0,
      #               0, 0, 0, 1, rhoattrit - .Machine$double.eps,
      #               0, 0, 0, rhoattrit - .Machine$double.eps, 1), nrow = 5, byrow = T))
      chol(matrix(c(sigeps00^2, 0, 0, 0, 0,
                    0, sigeps10^2, 0, 0, 0,
                    0, 0, sigeps11^2, 0, 0,
                    0, 0, 0, 1, rhoattrit + ifelse(rhoattrit > 0, -.Machine$double.eps, .Machine$double.eps),
                    0, 0, 0, rhoattrit + ifelse(rhoattrit > 0, -.Machine$double.eps, .Machine$double.eps), 1), nrow = 5, byrow = T))
    df <- data.frame(y00 = errors %*% matrix(c(1, 0, 0, beta1, beta2), nrow = 5),
                     y10 = errors %*% matrix(c(0, 1, 0, beta3, beta4), nrow = 5),
                     y11 = errors %*% matrix(c(0, 0, 1, beta3, beta4), nrow = 5) + tau,
                     phis0 = pnorm(errors[,4]), phis1 = pnorm(errors[,4]),
                     attrit10 = pnorm(errors[,4]) > 1 - pattrit0,
                     attrit11 = pnorm(errors[,5]) > 1 - pattrit1,
                     d = rep(c(0, 1), nsims * nobs / 2), i = rep(1:nobs, nsims),
                     sim = rep(1:nsims, each = nobs)) %>%
      mutate(y_0 = y00, y_1 = y10 * (1 - d) + y11 * d,
             y_1 = ifelse(d == 1 & attrit11, NA, ifelse(d == 0 & attrit10, NA, y_1)))
    
    df2 <- df %>% select(i, sim, d, y_0, y_1) %>%
      gather(t, y, -i, -sim, -d) %>% mutate(t = gsub("y_", "", t, fixed = T) %>% as.numeric) %>%
      mutate(attrit = is.na(y))
    df3 <- df2 %>% split(df2$sim)
    
    ols <- df3 %>%
      sapply(function(dfsim) lm(y ~ factor(d) * factor(t), data = dfsim)$coef["factor(d)1:factor(t)1"])
    # ols2 <- df3 %>%
    #   sapply(function(dfsim) felm(y ~ factor(d) * factor(t) | i + t, data = dfsim)$coef["factor(d)1:factor(t)1",1])
    balanced <- df2 %>% group_by(i, sim) %>% filter(all(!is.na(y))) %>% ungroup %>% split(.$sim) %>%
      sapply(function(dfsim) lm(y ~ factor(d) * factor(t), data = dfsim)$coef["factor(d)1:factor(t)1"])
    boundsdf <- df2 %>% group_by(i, sim) %>% mutate(dy = y - lag(y)) %>% ungroup %>% filter(t == 1) %>% split(.$sim)
    # manski <- boundsdf %>%
    #   sapply(function(dfsim) lm(y ~ factor(d) * factor(t), data = dfsim)$coef["factor(d)1:factor(t)1"])
    i <- 1
    lee <- boundsdf %>%
      sapply(function(dfsim) {
        # print(i); i <<- i + 1
        attrittreat <- sum(dfsim$attrit * dfsim$d)
        attritcontrol <- sum(dfsim$attrit * (1 - dfsim$d))
        if(attrittreat == attritcontrol) { dfsimlb <- dfsim; dfsimub <- dfsim }
        if(attrittreat > attritcontrol) {
          dfsimlb <- arrange(dfsim, d, dy)[-(1:(attrittreat - attritcontrol)),]
          dfsimub <- arrange(dfsim, d, -dy)[-(1:(attrittreat - attritcontrol)),]
        }
        if(attritcontrol > attrittreat) {
          dfsimlb <- arrange(dfsim, -d, -dy)[-(1:(attritcontrol - attrittreat)),]
          dfsimub <- arrange(dfsim, -d, dy)[-(1:(attritcontrol - attrittreat)),]
        }
        c(lm(dy ~ factor(d), data = dfsimlb)$coef["factor(d)1"],
          lm(dy ~ factor(d), data = dfsimub)$coef["factor(d)1"])
      })
    leelb <- lee[1,]
    leeub <- lee[2,]
    
    olsrange <- quantile(ols, c(.05, .95))
    olsest <- median(ols)
    balancedrange <- quantile(balanced, c(.05, .95))
    balancedest <- median(balanced)
    leerange <- c(quantile(leelb, .05), quantile(leeub, .95))
    leeest <- c(median(leelb), median(leeub))
    ggplot() +
      geom_hline(yintercept = tau, col = "red", linetype = 2) +
      annotate(geom = "label", x = -0.5, y = tau, col = "red", label = "True effect", size = 5) +
      annotate(geom = "errorbar", x = -1, ymin = olsrange[1], ymax = olsrange[2], width = 0.1) +
      annotate(geom = "point", x = -1, y = olsest) +
      annotate(geom = "errorbar", x = 0, ymin = balancedrange[1], ymax = balancedrange[2], width = 0.1) +
      annotate(geom = "point", x = 0, y = balancedest) +
      annotate(geom = "errorbar", x = 1, ymin = leerange[1], ymax = leerange[2], width = 0.1) +
      annotate(geom = "segment", x = 1, xend = 1, y = leeest[1], yend = leeest[2], size = 6) +
      theme_bw() + theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
                         axis.text.x = element_text(size = 18)) +
      scale_x_continuous(breaks = c(-1, 0, 1), labels = c("OLS", "Drop\nattriters", "Lee\nbounds"), name = NULL) +
      scale_y_continuous(name = NULL)
  })

  output$esplot <- renderPlot({
    genplotsim()
  })
}

shinyApp(ui, server)

