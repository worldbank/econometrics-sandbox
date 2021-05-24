library(dplyr)
library(ggplot2)
library(shiny)
library(shinydashboard)

# STEP 1: NON-LINEARITIES
# USE STEP FUNCTION (MPC0 --> MPC1) FOR SPECFICIATION OF f W/ KNOWN THRESHOLD
# REQUIRED: f(Size): T, MPC1, MPC2 [NOTE: THESE SHOULDN'T AFFECT SOLUTION, JUST GRAPH?]
# f(Size) = MPC1 + 1{Size > T} * (MPC2 - MPC1) * (Size - T) / Size
# OPTIMIZED: Size1*, Size2*, N0*, N1*, N2*
# REPORTED: Var*[beta2hat](Size1*, Size2*) / (sig2eps / N)
# PLOTTED:
# 1) f(Size) (always show 0, 3*T, Size2* + T)
# 2) Var*[beta2hat](Size1, Size2*) / (sig2eps / N)
# 3) Var*[beta2hat](Size1*, Size2) / (sig2eps / N)

# NOTE: INSTEAD JUST PROVIDE CALCULATOR FOR POWER FOR beta2 FOR f(size) = size^rho, ONLY OPTIMIZE N0/N1/N2

var_nl <- function(v) {
  stopifnot(all(c("f", "size_uct1", "size_uct2", "sdeps", "n0", "nuct1") %in% names(v)))
  stopifnot("Both UCT must have positive size" = v$size_uct2 > 0)
  if(v$n0 + v$nuct1 > 1 | v$n0 < 0 | v$nuct1 < 0) return(Inf)
  v$nuct2 <- 1 - (v$n0 + v$nuct1)
  var_beta <- (1/v$nuct2) * (1/v$size_uct2^2) + (1/v$nuct1) * (1/v$size_uct1^2) + (1/v$n0) * ((v$size_uct2 - v$size_uct1)^2 / (v$size_uct1^2 * v$size_uct2^2))
  var_beta <- v$sdeps^2 * var_beta * (1/(v$f(v$size_uct2) - v$f(v$size_uct1))^2)
  return(var_beta)
}

var_nl_reluct <- function(v) {
  stopifnot(all(c("f", "size_uct1", "size_uct2", "sdeps", "n0", "nuct1") %in% names(v)))
  stopifnot("Both UCT must have positive size" = v$size_uct2 > 0)
  if(v$n0 + v$nuct1 > 1 | v$n0 < 0 | v$nuct1 < 0) return(Inf)
  v$nuct2 <- 1 - (v$n0 + v$nuct1)
  var_beta <- (1/v$nuct2) * (1/v$size_uct2^2) + (1/v$nuct1) * (1/v$size_uct1^2) + (1/v$n0) * ((v$size_uct2 - v$size_uct1)^2 / (v$size_uct1^2 * v$size_uct2^2))
  var_beta <- v$sdeps^2 * var_beta * (v$size_uct2 / (v$size_uct2 - v$size_uct1))^2
  return(var_beta)
}

opt_nl <- function(v) {
  stopifnot(identical(names(v), c("f", "size_uct1", "size_uct2", "sdeps", "n0", "nuct1")))
  n <- 2*v$size_uct2
  v$n0 <- (v$size_uct2 - v$size_uct1) / n
  v$nuct1 <- v$size_uct2 / n
  return(v)
}

opt_nlnum <- function(v) {
  # v list OF VARAIBLES BELOW
  stopifnot(identical(names(v), c("f", "size_uct1", "size_uct2", "sdeps", "n0", "nuct1")))
  get_sol(v = v, default_v = c("n0" = 0.33, "nuct1" = 0.33),
          req_v = c("f", "size_uct1", "size_uct2", "sdeps"),
          f = var_nl)
}


# STEP 2: BENCHMARKING (SIMPLE CASE)
# REQUIRED: Size_UCT / Size_TUP, omega
# OPTIMIZED: N0*, N1*, N2*
# REPORTED: Var*[betahat_TUP-CT & betahat_TUP](omega; Size_UCT / Size_TUP), N0*, N1*, N2*
# PLOTTED:
# 1) Var*[betahat_TUP-CT & betahat_TUP](omega \in (0, 1); Size_UCT / Size_TUP) [TWO LINES]
# 2) N0* & N1* & N2* [betahat_TUP-CT & betahat_TUP](omega \in (0, 1); Size_UCT / Size_TUP) [THREE LINES]

var_bench <- function(v) {
  stopifnot(all(c("size_tup", "size_uct", "sdeps", "n0", "nuct") %in% names(v)))
  stopifnot(v$size_tup > 0 & v$size_uct > 0 & v$sdeps > 0)
  if(v$n0 + v$nuct > 1 | v$n0 < 0 | v$nuct < 0) return(c("var_betatup" = Inf, "var_betatupce" = Inf))
  v$ntup <- 1 - (v$n0 + v$nuct)
  var_betatup <- v$sdeps^2 * (1/v$size_tup^2) * (1/v$ntup + 1/v$n0)
  # FIX BELOW
  var_betatupce <- 1/v$ntup + ((1/v$nuct) * (v$size_tup^2/v$size_uct^2)) +
    ifelse(v$size_uct == v$size_tup, 0, ((1/v$n0) * (1 - (v$size_tup/v$size_uct))^2))
  var_betatupce <- v$sdeps^2 * (1/v$size_tup^2) * var_betatupce
  var_betauct <- v$sdeps^2 * (1/v$size_uct^2) * (1/v$nuct + 1/v$n0)
  return(c("var_betatup" = var_betatup, "var_betatupce" = var_betatupce, "var_betauct" = var_betauct))
}

var_benchomega <- function(v) {
  varv <- var_bench(v)
  v$omega3 <- 1 - v$omega1 - v$omega2
  if(isTRUE(all.equal(v$omega3, 0))) v$omega3 <- 0
  case_when(
    v$omega1 == 1 ~ varv["var_betatup"],
    v$omega2 == 1 ~ varv["var_betatupce"],
    v$omega3 == 1 ~ varv["var_betauct"],
    v$omega1 == 0 & v$omega2 > 0 & v$omega3 > 0 ~ v$omega2 * varv["var_betatupce"] + v$omega3 * varv["var_betauct"],
    v$omega1 > 0 & v$omega2 == 0 & v$omega3 > 0 ~ v$omega1 * varv["var_betatup"] + v$omega3 * varv["var_betauct"],
    v$omega1 > 0 & v$omega2 > 0 & v$omega3 == 0 ~ v$omega1 * varv["var_betatup"] + v$omega2 * varv["var_betatupce"],
    T ~ v$omega1 * varv["var_betatup"] + v$omega2 * varv["var_betatupce"] + v$omega3 * varv["var_betauct"]
  )
}

opt_bench <- function(v) {
  # v list OF VARAIBLES BELOW
  stopifnot(identical(names(v), c("omega1", "omega2", "size_tup", "size_uct", "sdeps", "n0", "nuct")))
  stopifnot("Weights must be non-negative and sum to less than or equal to 1." = v$omega1 >= 0 & v$omega2 >= 0 & v$omega1 + v$omega2 <= 1)
  n <- sqrt(v$omega1 + v$omega2 * (1 - (v$size_tup/v$size_uct))^2 + (1 - v$omega1 - v$omega2) * (v$size_tup/v$size_uct)^2) +
    sqrt(1 - v$omega1) * (v$size_tup/v$size_uct) + sqrt(v$omega1 + v$omega2)
  v$n0 <- sqrt(v$omega1 + v$omega2 * (1 - (v$size_tup/v$size_uct))^2 + (1 - v$omega1 - v$omega2) * (v$size_tup/v$size_uct)^2) / n
  v$nuct <- sqrt(1 - v$omega1) * (v$size_tup/v$size_uct) / n
  return(v)
}

opt_benchnum <- function(v) {
  # v list OF VARAIBLES BELOW
  stopifnot(identical(names(v), c("omega1", "omega2", "size_tup", "size_uct", "sdeps", "n0", "nuct")))
  get_sol(v = v, default_v = c("n0" = 0.33, "nuct" = 0.33),
          req_v = c("omega1", "omega2", "size_tup", "size_uct", "sdeps"),
          f = var_benchomega)
}


# STEP 3: BENCHMARKING (MZ18 CASE)
# REQUIRED: omega, Size_TUPmax / Size_TUPmin (greater than 1, less than Inf)
# OPTIMIZED: Size_UCT2 / Size_TUPmean, Size_UCT1 / Size_TUPmean, N0*, N1*, N2*
# REPORTED: Var*[betahat_TUP-CT & betahat_TUP](omega), N0*, N1*, N2*
# PLOTTED:
# 1) Var*[betahat_TUP-CT & betahat_TUP](omega \in (0, 1)) [TWO LINES]
# 2) N0* & N1* & N2* [betahat_TUP-CT & betahat_TUP](omega \in (0, 1)) [THREE LINES]

# INCREASE SIZE BUTTON, FUNCTIONALITY VARIES BY CHOICE OF rho...

var_mz18 <- function(v) {
  stopifnot(all(c("size_tupmin", "size_tupmax", "rho", "size_uct1", "size_uct2", "sdeps", "n0", "nuct1", "nuct2") %in% names(v)))
  # stopifnot(v$size_tupmin > 0 & v$size_tupmax > v$size_tupmin & v$rho != 0 & v$size_uct2 >= v$size_uct1)
  stopifnot(v$size_tupmin > 0 & v$size_tupmax > v$size_tupmin & v$rho != 0 & v$sdeps > 0)
  if(is.na(v$n0) | is.na(v$nuct1) | is.na(v$nuct2)) return(c("var_betatup" = NA, "var_betatupce" = NA, "var_betareluct" = NA))
  if(v$n0 + v$nuct1 + v$nuct2 > 1 | v$n0 < 0 | v$nuct1 < 0 | v$nuct2 < 0) return(c("var_betatup" = Inf, "var_betatupce" = Inf, "var_betareluct" = Inf))
  v$ntup <- 1 - (v$n0 + v$nuct1 + v$nuct2)
  var_betatup <- v$sdeps^2 * (4/(v$size_tupmin + v$size_tupmax)^2) * (1/v$ntup + 1/v$n0)
  getvar_betatupce <- function(s) {
    A31 <- s^v$rho - v$size_uct1^v$rho
    A21 <- v$size_uct2^v$rho - v$size_uct1^v$rho
    A23 <- v$size_uct2^v$rho - s^v$rho
    vardelta <- 1/v$ntup + ((1/v$nuct2) * (s / v$size_uct2)^2 * (A31 / A21)^2) + ((1/v$nuct1) * (s / v$size_uct1)^2 * (A23 / A21)^2) +
      ((1/v$n0) * (A21 - (s/v$size_uct2)*A31 - (s/v$size_uct1)*A23)^2 / A21^2)
    vardelta
  }
  var_betatupce <- v$sdeps^2 * (4/(v$size_tupmin + v$size_tupmax)^2) * max(getvar_betatupce(seq(v$size_tupmin, v$size_tupmax, length.out = 1e2)))
  var_betareluct <- (1/v$nuct2) * (1/v$size_uct2^2) + (1/v$nuct1) * (1/v$size_uct1^2) +
    (1/v$n0) * ((v$size_uct2 - v$size_uct1)^2 / (v$size_uct2^2 * v$size_uct1^2))
  var_betareluct <- v$sdeps^2 * var_betareluct * (max(v$size_uct1, v$size_uct2))^2 / (v$size_uct2 - v$size_uct1)^2
  return(c("var_betatup" = var_betatup, "var_betatupce" = var_betatupce, "var_betareluct" = var_betareluct))
}

var_mz18omega <- function(v) {
  varv <- var_mz18(v)
  v$omega3 <- 1 - v$omega1 - v$omega2
  if(isTRUE(all.equal(v$omega3, 0))) v$omega3 <- 0
  case_when(
    v$omega1 == 1 ~ varv["var_betatup"],
    v$omega2 == 1 ~ varv["var_betatupce"],
    v$omega3 == 1 ~ varv["var_betareluct"],
    v$omega1 == 0 & v$omega2 > 0 & v$omega3 > 0 ~ v$omega2 * varv["var_betatupce"] + v$omega3 * varv["var_betareluct"],
    v$omega1 > 0 & v$omega2 == 0 & v$omega3 > 0 ~ v$omega1 * varv["var_betatup"] + v$omega3 * varv["var_betareluct"],
    v$omega1 > 0 & v$omega2 > 0 & v$omega3 == 0 ~ v$omega1 * varv["var_betatup"] + v$omega2 * varv["var_betatupce"],
    T ~ v$omega1 * varv["var_betatup"] + v$omega2 * varv["var_betatupce"] + v$omega3 * varv["var_betareluct"]
  )
}

opt_mz18 <- function(v) {
  # v list OF VARAIBLES BELOW
  stopifnot(identical(names(v), c("size_tupmin", "size_tupmax", "omega1", "omega2", "rho", "size_uct1", "size_uct2", "sdeps", "n0", "nuct1", "nuct2")))
  # CHECK CONDITIONS ON sizetup_min AND size_tupmax (SHOULD HAPPEN EARLIER)
  v$omega3 <- 1 - v$omega1 - v$omega2
  if(isTRUE(all.equal(v$omega3, 0))) v$omega3 <- 0
  stopifnot(v$size_tupmin > 0 & v$size_tupmax > v$size_tupmin & v$sdeps > 0)
  stopifnot("Weights must be non-negative and sum to less than or equal to 1." = v$omega1 >= 0 & v$omega2 >= 0 & v$omega3 >= 0)
  
  get_sol(v = v, default_v = c("n0" = 0.25, "nuct1" = 0.25, "nuct2" = 0.25),
          req_v = c("size_uct1", "size_uct2", "size_tupmin", "size_tupmax", "omega1", "omega2", "omega3", "rho", "sdeps"),
          f = var_mz18omega)
}

# RETURNS MINIMIZERS OF f, HOLDING FIXED PROVIDED NON-NULL ELEMENTS OF v (REQUIREMENT THAT req_v NON-NULL), WITH INITIAL
# VALUES default_v FOR NULL ELEMENTS
# get_sol(list("x" = 0, "y" = 0, "a" = 1), default_v = c("x" = 0, "y" = 0), req_v = c("a"), f = function(v) { x <- v$x; y <- v$y; a <- v$a; (x + y + a)^2 + x^2 }) --> (0, 0, 1)
# get_sol(list("x" = NULL, "y" = 0, "a" = 1), default_v = c("x" = 0, "y" = 0), req_v = c("a"), f = function(v) { x <- v$x; y <- v$y; a <- v$a; (x + y + a)^2 + x^2 }) --> (-0.5, 0, 1)
# get_sol(list("x" = NULL, "y" = NULL, "a" = 1), default_v = c("x" = 0, "y" = 0), req_v = c("a"), f = function(v) { x <- v$x; y <- v$y; a <- v$a; (x + y + a)^2 + x^2 }) --> (0, -1, 1)
get_sol <- function(v, default_v, req_v, f) {
  # names(v) has no duplicates
  stopifnot(identical(names(v) %>% sort, unique(names(v)) %>% sort))
  # names(default_v) and req_v must union to names(v) w/ no overlap
  stopifnot(identical(c(names(default_v), req_v) %>% sort, names(v) %>% sort))
  
  opt_v <- names(v)[sapply(v, function(vi) {length(vi) == 0})]
  fix_v <- names(v)[sapply(v, function(vi) {length(vi) == 1})]
  
  stopifnot(all(names(v) %in% c(opt_v, fix_v)))
  stopifnot(all(req_v %in% fix_v))
  if(length(opt_v) == 0) {
    sol_v <- v
  } else {
    # CURRENTLY HANDLING ERROR FROM OPTIM WHEN INITIAL VALUES DON'T EVALUATE
    # COMES UP WHEN GRAPHING LIMIT CASES LIKE 0 SIZE SMALL UCT FOR RHO = -1 FOR MZ
    # SOLUTION: RETURN NAs IN THESE CASES (MEANS NEED TO HANDLE MISSINGS IN INPUTS TO var_mz18)
    sol_v <- tryCatch({
      optim(default_v[opt_v], function(theta) {
        names(theta) <- opt_v
        sols_v <- f(c(theta, unlist(v[fix_v]))[names(v)] %>% as.list)
        return(sols_v)
      }, control = list("maxit" = 5e3))$par
    },
    error = function(cond){ return(NULL)})
    
    if(length(sol_v) == 0) {
      sol_v <- c(setNames(rep(NA, length(opt_v)), opt_v), unlist(v[fix_v]))[names(v)] %>% as.list
    } else {
      names(sol_v) <- opt_v
      sol_v <- c(sol_v, unlist(v[fix_v]))[names(v)] %>% as.list
    }
  }
  return(sol_v)
}

ui <- fluidPage(
  titlePanel("Power Dashboard: Cash transfer size, nonlinearities, and benchmarking"),
  withMathJax(),
  tabsetPanel(
    tabPanel("Instructions",
             fluidRow(
               column(12,
                      helpText(paste0("This dashboard is a companion to Kondylis & Loeser (2021) and enables one to calculate the variance of estimators in three RCT designs ",
                                      "related to cash transfer size, nonlinearities, and benchmarking. The dashboard calculates the variance minimizing allocation ",
                                      "of subjects across a control group and treatment arms, and graphs how this variance changes as one changes key parameters. In each case,",
                                      "the dashboard enables one to vary the size of cash transfers used, and to calculate required sample sizes and minimum detectable effects ",
                                      "based on these variances. Variances ",
                                      "are calculated assuming homoskedasticity (constant variance of errors that does not change as parameters are changed), and variances ",
                                      "are reported relative to the variance of errors divided by the number of observations. The paper provides expressions for variance under ",
                                      "more general cases.")),
                      h4("Nonlinearities"),
                      helpText(paste0("The first case concerns testing whether the effects of cash transfers (UCT) are linear, ",
                                      "i.e., they are constant per unit of transfer with respect to transfer size. To do this, the researcher randomizes individuals across a control ",
                                      "group, a small UCT group, and a large UCT group. They then estimate whether the effects of the small and large UCT are equal per unit ",
                                      "of transfer. For evaluating statistical power in this case, it is important to specify the nonlinear model one is interested in ",
                                      "testing, and this dashboard reports the variance of the coefficient on the included nonlinear term. The relationship between the variance of ",
                                      "this coeffiicent and transfer sizes may vary with the nonlinear model, but the variance minimizing share of individuals allocated to each arm ",
                                      "does not depend on the nonlinear model selected with this setup. ",
                                      "Two nonlinear models are considered. The first is a power model, in which the effect of cash transfers per unit of transfer ",
                                      "equals a constant plus transfer size raised to a power. The second is a threshold model, which is consistent with a constant marginal ",
                                      "effect of cash transfers until a threshold T, above which the marginal effect discretely jumps. To facilitate interpretation, ",
                                      "minimum detectable effects are reported (and are specified) as the effect of the large UCT relative to the small UCT, scaled ",
                                      "per unit of transfer.")),
                      h4("Benchmarking"),
                      helpText(paste0("The second case concerns benchmarking a fixed program against UCT, i.e., estimating whether the effect of the program per unit ",
                                      "of cost is equal to the effect of UCT per unit of cost. To do this, the researcher randomizes individuals across a control group, a UCT, and ",
                                      "the program. They then estimate three effects. First, they estimate the effect of the program by comparing it to the control group ",
                                      "(\"Program Effect\"). Second, they estimate the effect of the program relative to the effect of the UCT, scaled by the relative cost of ",
                                      "the program (\"Program Effect Relative to Cash Benchmark\"). Third, they estimate the effect of the UCT by comparing it to the control group ",
                                      "(\"UCT Effect\"). The dashboard reports each of ",
                                      "these variances, and chooses the assignment of subjects across arms to minimize a weighted average of the variance of these three ",
                                      "estimates (given selected weights).")),
                      h4("Benchmarking with cost uncertainty"),
                      helpText(paste0("The third case concerns benchmarking a program with uncertain costs against UCT, while allowing for the possibility that ",
                                      "the effects of UCT per unit of cost may not be constant. To do this, the researcher randomizes individuals across a control group, ",
                                      "a small UCT, a large UCT, and the program. They then estimate three effects. First, they estimate the effect of the program by ",
                                      "comparing it to the control group. Second, following McIntosh & Zeitlin (2020), they estimate the effect of the program relative to ",
                                      "the predicted effect of a UCT of cost equal to the realized cost of the program. To do so, one must specify a non-constant model of ",
                                      "the effect of UCT; we consider a power model as in the nonlinear case, with the parametrization \\(\\rho = -1\\) corresponding to ",
                                      "McIntosh & Zeitlin (2020). Third, they estimate the effect of the additional transfers under the large UCT relative to the small UCT, ",
                                      "in both cases scaled per unit of transfer. The dashboard reports the first variance, the maximum of the second variance across the range of ",
                                      "possible costs of the program, and the third variance. For the second variance, one must therefore specify the range of anticipated costs. ",
                                      "Note that we implicitly assume that the effect of the program is unaffected by the realized cost. ",
                                      "The dashboard chooses the assignment of subjects to minimize a weighted average of the first variance, the ",
                                      "maximum of the second variance across potential cost realizations, and the third variance (given selected weights)."))
               )
             )
    ),
    tabPanel("Nonlinearities",
             titlePanel("Nonlinearities"),
             fluidRow(
               column(12,
                      helpText("\\[ Y_{i} = \\beta_{0} + \\beta_{1} \\text{Small UCT}_{i} + \\beta_{2} \\text{Large UCT}_{i} + \\epsilon_{i} \\]"),
                      helpText("\\[ Y_{i} = \\delta_{0} + \\delta_{1} \\text{Size}_{a(i)} + \\delta_{2} \\text{Size}_{a(i)} f(\\text{Size}_{a(i)}) + \\epsilon_{i} \\]"),
                      helpText(paste0("\\[ \\left( \\begin{array}{c} \\beta_{2}' \\\\ \\delta_{2} \\end{array} \\right) \\equiv ",
                                      "\\left( \\begin{array}{c} \\frac{\\beta_{2} - \\beta_{1}}{\\text{Size}_{\\text{Large UCT}} - \\text{Size}_{\\text{Small UCT}}}",
                                      " - \\frac{\\beta_{1}}{\\text{Size}_{\\text{Small UCT}}} \\\\ \\delta_{2} \\end{array} \\right) \\equiv ",
                                      "\\left( \\begin{array}{c} \\text{Large UCT marginal effect relative to Small UCT} \\\\ \\text{Nonlinear coefficient} \\end{array} \\right)",
                                      "\\]")),
                      helpText("\\[ \\min_{N_{\\text{Control}}, N_{\\text{Small UCT}}, N_{\\text{Large UCT}}} \\text{Var}[\\hat{\\delta}_{2}] \\]"),
                      helpText("\\[ N_{\\text{Control}} + N_{\\text{Small UCT}} + N_{\\text{Large UCT}} \\leq N \\]"),
                      selectInput("nl_fname", "Function \\( f \\)", c("Power", "Threshold")),
                      uiOutput("nl_fexpression"),
                      fluidRow(
                        column(6, textInput("nl_param1", uiOutput("nl_paramname"), 1)),
                        column(6, textInput("nl_sdeps", "\\( \\sqrt{\\text{Var}[\\epsilon]} \\): Std. dev. of error", 1))
                      )
               )
             ),
             fluidRow(
               column(6,
                      h4("Panel 1 Parameters"),
                      textInput("nl_panel1_size_uct1", "\\( \\text{Size}_{\\text{Small UCT}} \\): Size of Small UCT", 1),
                      textInput("nl_panel1_size_uct2", "\\( \\text{Size}_{\\text{Large UCT}} \\): Size of Large UCT", 2),
                      h4("Optimal"),
                      tableOutput("nl_panel1table")
               ),
               column(6,
                      h4("Panel 2 Parameters"),
                      textInput("nl_panel2_size_uct1", "\\( \\text{Size}_{\\text{Small UCT}} \\): Size of Small UCT", 1),
                      textInput("nl_panel2_size_uct2", "\\( \\text{Size}_{\\text{Large UCT}} \\): Size of Large UCT", 4),
                      h4("Optimal"),
                      tableOutput("nl_panel2table")
               )
             ),
             fluidRow(
               column(12,
                      tabsetPanel(
                        tabPanel("Power calculation (Panel 1 Parameters)",
                                 fluidRow(
                                   column(6, textInput("nl_power", "Power", 0.8)),
                                   column(6, textInput("nl_alpha", "Alpha", 0.05))
                                 ),
                                 fluidRow(
                                   column(6,
                                          h4("Minimum Detectable Effect (MDE)"),
                                          textInput("nl_samplesize", "Sample Size \\( N \\)", 1000),
                                          tableOutput("nl_panel1_table_mde")
                                   ),
                                   column(6,
                                          h4("Required Sample Size"),
                                          textInput("nl_mde_betareluct", "MDE (Large UCT marg. effect rel. to Small UCT)", 0.1),
                                          uiOutput("nl_delta"),
                                          tableOutput("nl_panel1_table_samplesize")
                                   )
                                 )
                        ),
                        tabPanel("# of obs. by \\( \\text{Size}_{\\text{Small UCT}} \\)",
                                 fluidRow(
                                   column(6, plotOutput("nl_panel1_plot_n_size_uct1")),
                                   column(6, plotOutput("nl_panel2_plot_n_size_uct1"))
                                 )
                        ),
                        tabPanel("# of obs. by \\( \\text{Size}_{\\text{Large UCT}} \\)",
                                 fluidRow(
                                   column(6, plotOutput("nl_panel1_plot_n_size_uct2")),
                                   column(6, plotOutput("nl_panel2_plot_n_size_uct2"))
                                 )
                        ),
                        tabPanel("Variance by \\( \\text{Size}_{\\text{Small UCT}} \\)",
                                 fluidRow(
                                   column(6, plotOutput("nl_panel1_plot_var_size_uct1")),
                                   column(6, plotOutput("nl_panel2_plot_var_size_uct1"))
                                 )
                        ),
                        tabPanel("Variance by \\( \\text{Size}_{\\text{Large UCT}} \\)",
                                 fluidRow(
                                   column(6, plotOutput("nl_panel1_plot_var_size_uct2")),
                                   column(6, plotOutput("nl_panel2_plot_var_size_uct2"))
                                 )
                        )
                      )
               )
             )
    ),
    tabPanel("Benchmarking",
             titlePanel("Benchmarking"),
             fluidRow(
               column(12,
                      helpText("\\[ Y_{i} = \\beta_{0} + \\beta_{1} \\text{UCT}_{i} + \\beta_{2} \\text{Program}_{i} + \\epsilon_{i} \\]"),
                      helpText(paste0("\\[ \\left( \\begin{array}{c} \\beta_{2}' \\\\ \\delta_{2}' \\\\ \\beta_{1}'  \\end{array} \\right) \\equiv ",
                                      "\\left( \\begin{array}{c} \\frac{1}{\\text{Size}_{\\text{Program}}} \\beta_{2} \\\\ ",
                                      "\\frac{1}{\\text{Size}_{\\text{Program}}} \\beta_{2} - \\frac{1}{\\text{Size}_{\\text{UCT}}} \\beta_{1} \\\\ ",
                                      "\\frac{1}{\\text{Size}_{\\text{UCT}}} \\beta_{1} ",
                                      "\\end{array} \\right) \\equiv \\left( \\begin{array}{c} \\text{Program effect} \\\\ ",
                                      "\\text{Program effect relative to cash benchmark} \\\\ \\text{UCT effect} \\end{array} \\right) \\]")),
                      helpText(paste0("\\[ \\min_{N_{\\text{Control}}, N_{\\text{UCT}}, N_{\\text{Program}}} \\omega_{1} \\text{Var}[\\hat{\\beta'_{2}}] + ",
                                      "\\omega_{2} \\text{Var}[\\hat{\\delta'_{2}}] + (1 - \\omega_{1} - \\omega_{2}) \\text{Var}[\\hat{\\beta'_{1}}] \\]")),
                      helpText("\\[ N_{\\text{Control}} + N_{\\text{UCT}} + N_{\\text{Program}} \\leq N \\]")
               )
             ),
             fluidRow(
               column(6, textInput("bm_size_tup", "\\( \\text{Size}_{\\text{Program}} \\): Anticipated cost", 1)),
               column(6, textInput("bm_sdeps", "\\( \\sqrt{\\text{Var}[\\epsilon]} \\): Std. dev. of error", 1))
             ),
             fluidRow(
               column(6,
                      h4("Panel 1 Parameters"),
                      textInput("bm_panel1_size_uct", "\\( \\text{Size}_{\\text{UCT}} \\): Size of UCT", 1),
                      textInput("bm_panel1_omega1", "\\( \\omega_{1} \\): Weight on variance of \\( \\hat{\\beta'_{2}} \\)", 0.5),
                      textInput("bm_panel1_omega2", "\\( \\omega_{2} \\): Weight on variance of \\( \\hat{\\delta'_{2}} \\)", 0.5),
                      h4("Optimal"),
                      tableOutput("bm_panel1table")
               ),
               column(6,
                      h4("Panel 2 Parameters"),
                      textInput("bm_panel2_size_uct", "\\( \\text{Size}_{\\text{UCT}} \\): Size of UCT", 2),
                      textInput("bm_panel2_omega1", "\\( \\omega_{1} \\): Weight on variance of \\( \\hat{\\beta'_{2}} \\)", 0.5),
                      textInput("bm_panel2_omega2", "\\( \\omega_{2} \\): Weight on variance of \\( \\hat{\\delta'_{2}} \\)", 0.5),
                      h4("Optimal"),
                      tableOutput("bm_panel2table")
               )
             ),
             fluidRow(
               column(12,
                      tabsetPanel(
                        tabPanel("Power calculation (Panel 1 Parameters)",
                                 fluidRow(
                                   column(6, textInput("bm_power", "Power", 0.8)),
                                   column(6, textInput("bm_alpha", "Alpha", 0.05))
                                 ),
                                 fluidRow(
                                   column(6,
                                          h4("Minimum Detectable Effect (MDE)"),
                                          textInput("bm_samplesize", "Sample Size \\( N \\)", 1000),
                                          tableOutput("bm_panel1_table_mde")
                                   ),
                                   column(6,
                                          h4("Required Sample Size"),
                                          textInput("bm_mde_betatup", "MDE (Program effect)", 0.5),
                                          textInput("bm_mde_betatupce", "MDE (Program effect rel. to cash benchmark)", 0.3),
                                          textInput("bm_mde_betauct", "MDE (UCT effect)", 0.5),
                                          tableOutput("bm_panel1_table_samplesize")
                                   )
                                 )
                        ),
                        tabPanel("# of obs. by \\( \\text{Size}_{\\text{UCT}} \\)",
                                 fluidRow(
                                   column(6, plotOutput("bm_panel1_plot_n_size_uct")),
                                   column(6, plotOutput("bm_panel2_plot_n_size_uct"))
                                 )
                        ),
                        tabPanel("# of obs. by \\( \\omega_{1} \\)",
                                 fluidRow(
                                   column(6, plotOutput("bm_panel1_plot_n_omega1")),
                                   column(6, plotOutput("bm_panel2_plot_n_omega1"))
                                 )
                        ),
                        tabPanel("# of obs. by \\( \\omega_{2} \\)",
                                 fluidRow(
                                   column(6, plotOutput("bm_panel1_plot_n_omega2")),
                                   column(6, plotOutput("bm_panel2_plot_n_omega2"))
                                 )
                        ),
                        tabPanel("Variance by \\( \\text{Size}_{\\text{UCT}} \\)",
                                 fluidRow(
                                   column(6, plotOutput("bm_panel1_plot_var_size_uct")),
                                   column(6, plotOutput("bm_panel2_plot_var_size_uct"))
                                 )
                        ),
                        tabPanel("Variance by \\( \\omega_{1} \\)",
                                 fluidRow(
                                   column(6, plotOutput("bm_panel1_plot_var_omega1")),
                                   column(6, plotOutput("bm_panel2_plot_var_omega1"))
                                 )
                        ),
                        tabPanel("Variance by \\( \\omega_{2} \\)",
                                 fluidRow(
                                   column(6, plotOutput("bm_panel1_plot_var_omega2")),
                                   column(6, plotOutput("bm_panel2_plot_var_omega2"))
                                 )
                        )
                      )
               )
             )
    ),
    tabPanel("Benchmarking with cost uncertainty",
             titlePanel("Benchmarking with cost uncertainty"),
             fluidRow(
               column(12,
                      helpText(paste0("\\[ Y_{i} = \\beta_{0} + \\beta_{1} \\text{Small UCT}_{i} + \\beta_{2} \\text{Large UCT}_{i} + ",
                                      "\\beta_{3} \\text{Program}_{i} + \\epsilon_{i} \\]")),
                      helpText(paste0("\\[ Y_{i} = \\delta_{0} + \\delta_{1} \\text{Size}_{a(i)} + \\delta_{2} \\text{Size}_{a(i)}^{1+\\rho} + ",
                                      "\\delta_{3} \\text{Program}_{i} + \\epsilon_{i} \\]")),
                      helpText(paste0("\\[ \\left( \\begin{array}{c} \\beta_{3}' \\\\ \\delta_{3}' \\\\ \\beta_{2}' - \\beta_{1}' \\end{array} \\right) \\equiv ",
                                      "\\left( \\begin{array}{c} \\frac{1}{(\\overline{\\text{Size}_{\\text{Program}}} + \\underline{\\text{Size}_{\\text{Program}}}) / 2} \\beta_{3} \\\\ ",
                                      "\\frac{1}{(\\overline{\\text{Size}_{\\text{Program}}} + \\underline{\\text{Size}_{\\text{Program}}}) / 2} \\delta_{3} \\\\ ",
                                      "\\frac{\\beta_{2} - \\beta_{1}}{\\text{Size}_{\\text{Large UCT}} - \\text{Size}_{\\text{Small UCT}}} - ",
                                      "\\frac{\\beta_{1}}{\\text{Size}_{\\text{Small UCT}}} ",
                                      "\\end{array} \\right) \\equiv \\left( \\begin{array}{c} \\text{Program effect} \\\\ ",
                                      "\\text{Program effect relative to cash benchmark} \\\\ ",
                                      "\\text{Large UCT marginal effect relative to Small UCT} \\end{array} \\right) \\]")),
                      # helpText(paste0("\\[ \\min_{N_{\\text{Control}}, N_{\\text{Small UCT}}, N_{\\text{Large UCT}}, N_{\\text{Program}}} \\left\\{ ",
                      #                 "\\omega_{1} \\frac{\\text{Var}[\\hat{\\beta}_{3}]}{\\mathbf{E}[\\text{Size}_{\\text{Program}}]^2} + ",
                      #                 "\\omega_{2} \\frac{1}{\\mathbf{E}[\\text{Size}_{\\text{Program}}]^2} \\max_{\\text{Size}_{\\text{Program}} \\in ",
                      #                 "[ \\underline{\\text{Size}_{\\text{Program}}}, \\overline{\\text{Size}_{\\text{Program}}} ]} ",
                      #                 "\\text{Var}[\\hat{\\delta}_{3}] \\right\\} \\]")),
                      helpText(paste0("\\[ \\min_{N_{\\text{Control}}, N_{\\text{Small UCT}}, N_{\\text{Large UCT}}, N_{\\text{Program}}} \\left\\{ ",
                                      "\\omega_{1} \\text{Var}[\\hat{\\beta_{3}'}] + ",
                                      "\\omega_{2} \\max_{\\text{Size}_{\\text{Program}} \\in ",
                                      "[ \\underline{\\text{Size}_{\\text{Program}}}, \\overline{\\text{Size}_{\\text{Program}}} ]} ",
                                      "\\text{Var}[\\hat{\\delta_{3}'}] + (1 - \\omega_{1} - \\omega_{2}) \\text{Var}[\\hat{\\beta_{2}'} - \\hat{\\beta_{1}'}] \\right\\} \\]")),
                      helpText(paste0("\\[ N_{\\text{Control}} + N_{\\text{Small UCT}} + N_{\\text{Large UCT}} + N_{\\text{Program}} \\leq N \\]")),
                      fluidRow(
                        column(6, textInput("mz_param1", "Parameter \\( \\rho \\)", -1)),
                        column(6, textInput("mz_sdeps", "\\( \\sqrt{\\text{Var}[\\epsilon]} \\): Std. dev. of error", 1))
                      ),
                      fluidRow(
                        column(6, textInput("mz_size_tupmin", "\\( \\underline{\\text{Size}_{\\text{Program}}} \\): Min. anticipated cost", 0.8)),
                        column(6, textInput("mz_size_tupmax", "\\( \\overline{\\text{Size}_{\\text{Program}}} \\): Max. anticipated cost", 1.2))
                      )
               )
             ),
             fluidRow(
               column(6,
                      h4("Panel 1 Parameters"),
                      textInput("mz_panel1_size_uct1", "\\( \\text{Size}_{\\text{Small UCT}} \\): Size of Small UCT", 0.8),
                      textInput("mz_panel1_size_uct2", "\\( \\text{Size}_{\\text{Large UCT}} \\): Size of Large UCT", 1.2),
                      textInput("mz_panel1_omega1", "\\( \\omega_{1} \\): Weight on variance of \\( \\hat{\\beta'_{3}} \\)", 0.5),
                      textInput("mz_panel1_omega2", "\\( \\omega_{2} \\): Weight on variance of \\( \\hat{\\delta'_{3}} \\)", 0.5),
                      h4("Optimal"),
                      tableOutput("mz_panel1table")
               ),
               column(6,
                      h4("Panel 2 Parameters"),
                      textInput("mz_panel2_size_uct1", "\\( \\text{Size}_{\\text{Small UCT}} \\): Size of Small UCT", 0.4),
                      textInput("mz_panel2_size_uct2", "\\( \\text{Size}_{\\text{Large UCT}} \\): Size of Large UCT", 1.6),
                      textInput("mz_panel2_omega1", "\\( \\omega_{1} \\): Weight on variance of \\( \\hat{\\beta'_{3}} \\)", 0.5),
                      textInput("mz_panel2_omega2", "\\( \\omega_{2} \\): Weight on variance of \\( \\hat{\\delta'_{3}} \\)", 0.5),
                      h4("Optimal"),
                      tableOutput("mz_panel2table")
               )
             ),
             fluidRow(
               column(12,
                      tabsetPanel(
                        tabPanel("Power calculation (Panel 1 Parameters)",
                                 fluidRow(
                                   column(6, textInput("mz_power", "Power", 0.8)),
                                   column(6, textInput("mz_alpha", "Alpha", 0.05))
                                 ),
                                 fluidRow(
                                   column(6,
                                          h4("Minimum Detectable Effect (MDE)"),
                                          textInput("mz_samplesize", "Sample Size \\( N \\)", 1000),
                                          tableOutput("mz_panel1_table_mde")
                                   ),
                                   column(6,
                                          h4("Required Sample Size"),
                                          textInput("mz_mde_betatup", "MDE (Program effect)", 0.5),
                                          textInput("mz_mde_betatupce", "MDE (Program effect rel. to cash benchmark)", 0.3),
                                          textInput("mz_mde_betareluct", "MDE (Large UCT marg. effect rel. to Small UCT)", 0.1),
                                          tableOutput("mz_panel1_table_samplesize")
                                   )
                                 )
                        ),
                        tabPanel("# of obs. by \\( \\text{Size}_{\\text{Small UCT}} \\)",
                                 fluidRow(
                                   column(6, plotOutput("mz_panel1_plot_n_size_uct1")),
                                   column(6, plotOutput("mz_panel2_plot_n_size_uct1"))
                                 )
                        ),
                        tabPanel("# of obs. by \\( \\text{Size}_{\\text{Large UCT}} \\)",
                                 fluidRow(
                                   column(6, plotOutput("mz_panel1_plot_n_size_uct2")),
                                   column(6, plotOutput("mz_panel2_plot_n_size_uct2"))
                                 )
                        ),
                        tabPanel("# of obs. by \\( \\omega_{1} \\)",
                                 fluidRow(
                                   column(6, plotOutput("mz_panel1_plot_n_omega1")),
                                   column(6, plotOutput("mz_panel2_plot_n_omega1"))
                                 )
                        ),
                        tabPanel("# of obs. by \\( \\omega_{2} \\)",
                                 fluidRow(
                                   column(6, plotOutput("mz_panel1_plot_n_omega2")),
                                   column(6, plotOutput("mz_panel2_plot_n_omega2"))
                                 )
                        ),
                        tabPanel("Variance by \\( \\text{Size}_{\\text{Small UCT}} \\)",
                                 fluidRow(
                                   column(6, plotOutput("mz_panel1_plot_var_size_uct1")),
                                   column(6, plotOutput("mz_panel2_plot_var_size_uct1"))
                                 )
                        ),
                        tabPanel("Variance by \\( \\text{Size}_{\\text{Large UCT}} \\)",
                                 fluidRow(
                                   column(6, plotOutput("mz_panel1_plot_var_size_uct2")),
                                   column(6, plotOutput("mz_panel2_plot_var_size_uct2"))
                                 )
                        ),
                        tabPanel("Variance by \\( \\omega_{1} \\)",
                                 fluidRow(
                                   column(6, plotOutput("mz_panel1_plot_var_omega1")),
                                   column(6, plotOutput("mz_panel2_plot_var_omega1"))
                                 )
                        ),
                        tabPanel("Variance by \\( \\omega_{2} \\)",
                                 fluidRow(
                                   column(6, plotOutput("mz_panel1_plot_var_omega2")),
                                   column(6, plotOutput("mz_panel2_plot_var_omega2"))
                                 )
                        )
                      )
               )
             )
    )
  )
)

server <- function(input, output) {
  nl_fname <- reactive(input$nl_fname)
  output$nl_fexpression <- renderUI({
    if(nl_fname() == "Power") nl_fout <- "\\( f(s) = s^{\\rho} \\)"
    if(nl_fname() == "Threshold") nl_fout <- "\\( f(s) = \\mathbf{1} \\{ s > T \\} \\frac{s - T}{s} \\)"
    nl_fout <- paste0("Selected function: ", nl_fout)
    withMathJax(
      tags$p(nl_fout)
    )
  })
  output$nl_paramname <- renderUI({
    if(nl_fname() == "Power") nl_paramnameout <- "\\( \\rho \\)"
    if(nl_fname() == "Threshold") nl_paramnameout <- "\\( T \\)"
    nl_paramnameout <- paste0("Parameter ", nl_paramnameout)
    withMathJax(
      tags$p(nl_paramnameout)
    )
  })
  
  nl_param1 <- reactive(as.numeric(input$nl_param1))
  nl_mde_betareluct <- reactive(as.numeric(input$nl_mde_betareluct))
  nl_size_uct1 <- reactive(as.numeric(input$nl_panel1_size_uct1))
  nl_size_uct2 <- reactive(as.numeric(input$nl_panel1_size_uct2))
  output$nl_delta <- renderUI({
    stopifnot("Large UCT must be larger than Small UCT" = nl_size_uct2() > nl_size_uct1())
    stopifnot("Both UCT must have positive size" = nl_size_uct2() > 0)
    nl_deltamde <- nl_mde_betareluct() * ((nl_size_uct2() - nl_size_uct1()) / nl_size_uct2()) * case_when(
      nl_fname() == "Power" ~  1 / (nl_size_uct2()^nl_param1() - nl_size_uct1()^nl_param1()),
      nl_fname() == "Threshold" ~ 1 / ((nl_size_uct2() >= nl_param1()) * ((nl_size_uct2() - nl_param1()) / nl_size_uct2()) +
                                          (nl_size_uct1() >= nl_param1()) * ((nl_size_uct1() - nl_param1()) / nl_size_uct1()))
    )
    # nl_deltaout <- paste0("Implied MDE (\\( \\delta_{2} \\)) = ", nl_deltamde)
    nl_deltaout <- paste0("Implied MDE (Nonlinear coefficient) = ", nl_deltamde)
    withMathJax(
      tags$p(nl_deltaout)
    )
  })
  
  nl_getfunction <- function(nl_param1, nl_fname) {
    f <- NULL
    if(nl_fname == "Power") {
      f <- function(x) x^nl_param1
    } else if(nl_fname == "Threshold") {
      f <- function(x) (x >= nl_param1) * (x - nl_param1) / x
    } else {
      stop("FUNCTION NOT PROPERLY SPECIFIED")
    }
    return(f)
  }
  
  nl_genpaneltable <- function(isize_uct1, isize_uct2) {
    f <- nl_getfunction(as.numeric(input$nl_param1), input$nl_fname)
    size_uct1 <- as.numeric(isize_uct1)
    size_uct2 <- as.numeric(isize_uct2)
    sdeps <- as.numeric(input$nl_sdeps)
    stopifnot("Large UCT must be larger than Small UCT" = size_uct2 > size_uct1)
    soln <- opt_nl(list("f" = f, "size_uct1" = size_uct1, "size_uct2" = size_uct2, "sdeps" = sdeps, "n0" = NULL, "nuct1" = NULL))
    tab <- data.frame(c("Fraction Control" = soln$n0, "Fraction Small UCT" = soln$nuct1, "Fraction Large UCT" = 1 - (soln$n0 + soln$nuct1),
                        # "Variance \\( \\hat{\\beta}_{2} \\)" = var_nl(soln))) %>%
                        "Var[Nonlinear Coefficient]" = var_nl(soln))) %>%
      as.matrix
    colnames(tab) <- " "
    tab
  }
  
  nl_genpowertable <- function(powervar, isize_uct1, isize_uct2) {
    f <- nl_getfunction(as.numeric(input$nl_param1), input$nl_fname)
    size_uct1 <- as.numeric(isize_uct1)
    size_uct2 <- as.numeric(isize_uct2)
    sdeps <- as.numeric(input$nl_sdeps)
    soln <- opt_nl(list("f" = f, "size_uct1" = size_uct1, "size_uct2" = size_uct2, "sdeps" = sdeps, "n0" = NULL, "nuct1" = NULL))
    varsoln <- var_nl_reluct(soln)
    alpha <- as.numeric(input$nl_alpha)
    power <- as.numeric(input$nl_power)
    N <- as.numeric(input$nl_samplesize)
    mde <- as.numeric(input$nl_mde_betareluct)
    if(powervar == "mde") {
      out <- unname((qnorm(1 - alpha/2) + qnorm(power)) * sqrt(varsoln / N))
    } else if(powervar == "n") {
      out <- unname((qnorm(1 - alpha/2) + qnorm(power))^2 * varsoln / (mde^2))
    }
    outnames <- c("Large UCT Marg. Effect Rel. to Small UCT")
    tab <- data.frame(setNames(out, outnames)) %>% as.matrix
    colnames(tab) <- " "
    tab
  }
  
  nl_genplot <- function(j, plotxvar, plotyvar, isize_uct1, isize_uct2, isize_uct1_other, isize_uct2_other) {
    f <- nl_getfunction(as.numeric(input$nl_param1), input$nl_fname)
    size_uct1 <- as.numeric(isize_uct1)
    size_uct2 <- as.numeric(isize_uct2)
    size_uct1_other <- as.numeric(isize_uct1_other)
    size_uct2_other <- as.numeric(isize_uct2_other)
    sdeps <- as.numeric(input$nl_sdeps)
    soln <- opt_nl(list("f" = f, "size_uct1" = size_uct1, "size_uct2" = size_uct2, "sdeps" = sdeps, "n0" = NULL, "nuct1" = NULL))
    varsoln <- var_nl(soln)
    soln_other <- opt_nl(list("f" = f, "size_uct1" = size_uct1_other, "size_uct2" = size_uct2_other, "sdeps" = sdeps, "n0" = NULL, "nuct1" = NULL))
    varsoln_other <- var_nl(soln_other)
    # GENERATE PLOT DATA
    if(plotxvar == "size_uct1") {
      v <- lapply(seq(0, size_uct2, length.out = 41), function(x) {
        soln_size_uct1 <- opt_nl(list("f" = f, "size_uct1" = x, "size_uct2" = size_uct2, "sdeps" = sdeps, "n0" = NULL, "nuct1" = NULL))
        varsoln_size_uct1 <- var_nl(soln_size_uct1)
        list("x" = x, "soln" = soln_size_uct1, "varsoln" = varsoln_size_uct1)
      })
    } else if(plotxvar == "size_uct2") {
      v <- lapply(seq(size_uct1, 2*max(2, size_uct2, size_uct2_other), length.out = 41), function(x) {
        soln_size_uct2 <- opt_nl(list("f" = f, "size_uct1" = size_uct1, "size_uct2" = x, "sdeps" = sdeps, "n0" = NULL, "nuct1" = NULL))
        varsoln_size_uct2 <- var_nl(soln_size_uct2)
        list("x" = x, "soln" = soln_size_uct2, "varsoln" = varsoln_size_uct2)
      })
    }
    if(plotyvar == "n") {
      df <- lapply(v, function(vi) {
        data.frame(x = vi$x, group = c("n0", "n1", "n2"),
                   y = c(vi[["soln"]]$n0, vi[["soln"]]$nuct1, 1 - (vi[["soln"]]$n0 +  vi[["soln"]]$nuct1)))
      }) %>% bind_rows()
    } else if(plotyvar == "var") {
      df <- lapply(v, function(vi) {
        data.frame(x = vi$x, group = c("var0"),
                   y = c(vi[["varsoln"]]))
      }) %>% bind_rows()
    }
    if(plotyvar == "n") {
      labels <- c(expression(N[Control]/N), expression(N[Small~UCT]/N), expression(N[Large~UCT]/N))
    } else if(plotyvar == "var") {
      labels <- c(expression(paste("Var[", hat("\u03B2")[2], "]")))
    }
    # GENERATE PLOT
    g <- ggplot(df, aes(x = x, y = y, group = group, col = group)) +
      geom_line() +
      geom_vline(xintercept = case_when(plotxvar == "size_uct1" ~ size_uct1, plotxvar == "size_uct2" ~ size_uct2), linetype = 2) +
      xlim(0, case_when(plotxvar == "size_uct1" ~ 2*max(size_uct2, size_uct2_other, 2),
                        plotxvar == "size_uct2" ~ 2*max(size_uct2, size_uct2_other, 2))) +
      # ylim(0, case_when(plotyvar == "n" ~ 1, plotyvar == "var" ~ 2 * max(varsoln, varsoln_other))) +
      xlab(case_when(plotxvar == "size_uct1" ~ expression(Size[Small~UCT]), plotxvar == "size_uct2" ~ expression(Size[Large~UCT]))) +
      # ylab(case_when(plotyvar == "n" ~ "Fraction of obs.", plotyvar == "var" ~ "Variance")) +
      scale_color_discrete(labels = labels) +
      theme_bw(base_size = 16) +
      theme(legend.title = element_blank(), legend.position = c(0.7, 0.8)) +
      ggtitle(paste0("Panel ", j))
    if(plotyvar == "n") {
      g <- g + ylim(0, 1) + ylab("Fraction of obs.")
    } else if(plotyvar == "var") {
      g <- g + scale_y_log10(limits = c(0.8 * min(df$y), 2 * max(varsoln, varsoln_other)), name = "Variance (log scale)")
    }
    return(g)
  }
  
  bm_genpaneltable <- function(isize_uct, iomega1, iomega2) {
    size_uct <- as.numeric(isize_uct)
    omega1 <- as.numeric(iomega1)
    omega2 <- as.numeric(iomega2)
    size_tup <- as.numeric(input$bm_size_tup)
    sdeps <- as.numeric(input$bm_sdeps)
    soln <- opt_bench(list("omega1" = omega1, "omega2" = omega2, "size_tup" = size_tup, "size_uct" = size_uct, "sdeps" = sdeps, "n0" = NULL, "nuct" = NULL))
    varsoln <- var_bench(soln)
    tab <- data.frame(c("Fraction Control" = soln$n0, "Fraction UCT" = soln$nuct, "Fraction Program" = 1 - (soln$n0 + soln$nuct),
                        "Var[Est. Program Effect]" = varsoln["var_betatup"] %>% unname,
                        "Var[Est. Program Effect Rel. to Benchmark]" = varsoln["var_betatupce"] %>% unname,
                        "Var[Est. UCT Effect]" = varsoln["var_betauct"] %>% unname)) %>%
      as.matrix
    colnames(tab) <- " "
    tab
  }
  
  bm_genpowertable <- function(powervar, isize_uct, iomega1, iomega2) {
    size_uct <- as.numeric(isize_uct)
    omega1 <- as.numeric(iomega1)
    omega2 <- as.numeric(iomega2)
    size_tup <- as.numeric(input$bm_size_tup)
    sdeps <- as.numeric(input$bm_sdeps)
    soln <- opt_bench(list("omega1" = omega1, "omega2" = omega2, "size_tup" = size_tup, "size_uct" = size_uct, "sdeps" = sdeps, "n0" = NULL, "nuct" = NULL))
    varsoln <- var_bench(soln)
    alpha <- as.numeric(input$bm_alpha)
    power <- as.numeric(input$bm_power)
    N <- as.numeric(input$bm_samplesize)
    mde <- c(as.numeric(input$bm_mde_betatup), as.numeric(input$bm_mde_betatupce), as.numeric(input$bm_mde_betauct))
    if(powervar == "mde") {
      out <- unname((qnorm(1 - alpha/2) + qnorm(power)) * sqrt(varsoln / N))
    } else if(powervar == "n") {
      out <- unname((qnorm(1 - alpha/2) + qnorm(power))^2 * varsoln / (mde^2))
    }
    outnames <- c("Program Effect", "Program Effect Rel. to Cash Benchmark", "UCT Effect")
    tab <- data.frame(setNames(out, outnames)) %>% as.matrix
    colnames(tab) <- " "
    tab
  }
  
  bm_genplot <- function(j, plotxvar, plotyvar, isize_uct, iomega1, iomega2, isize_uct_other, iomega1_other, iomega2_other) {
    size_uct <- as.numeric(isize_uct)
    omega1 <- as.numeric(iomega1)
    omega2 <- as.numeric(iomega2)
    size_uct_other <- as.numeric(isize_uct_other)
    omega1_other <- as.numeric(iomega1_other)
    omega2_other <- as.numeric(iomega2_other)
    size_tup <- as.numeric(input$bm_size_tup)
    sdeps <- as.numeric(input$bm_sdeps)
    soln <- opt_bench(list("omega1" = omega1, "omega2" = omega2, "size_tup" = size_tup, "size_uct" = size_uct, "sdeps" = sdeps, "n0" = NULL, "nuct" = NULL))
    varsoln <- var_bench(soln)
    soln_other <- opt_bench(list("omega1" = omega1_other, "omega2" = omega2_other, "size_tup" = size_tup, "size_uct" = size_uct_other, "sdeps" = sdeps, "n0" = NULL, "nuct" = NULL))
    varsoln_other <- var_bench(soln_other)
    # GENERATE PLOT DATA
    if(plotxvar == "omega1") {
      v <- lapply(seq(0.0125, 0.9875, 0.025), function(x) {
        soln_omega1 <- opt_bench(list("omega1" = x, "omega2" = (1 - x) * omega2 / (1 - omega1), "size_tup" = size_tup, "size_uct" = size_uct, "sdeps" = sdeps, "n0" = NULL, "nuct" = NULL))
        varsoln_omega1 <- var_bench(soln_omega1)
        list("x" = x, "soln" = soln_omega1, "varsoln" = varsoln_omega1)
      })
    } else if(plotxvar == "omega2") {
      v <- lapply(seq(0.0125, 0.9875, 0.025), function(x) {
        soln_omega2 <- opt_bench(list("omega1" = (1 - x) * omega1 / (1 - omega2), "omega2" = x, "size_tup" = size_tup, "size_uct" = size_uct, "sdeps" = sdeps, "n0" = NULL, "nuct" = NULL))
        varsoln_omega2 <- var_bench(soln_omega2)
        list("x" = x, "soln" = soln_omega2, "varsoln" = varsoln_omega2)
      })
    } else if(plotxvar == "size_uct") {
      v <- lapply(seq(min(0.5, size_uct, size_uct_other), 2*max(2, size_uct, size_uct_other), length.out = 41), function(x) {
        soln_size_uct <- opt_bench(list("omega1" = omega1, "omega2" = omega2, "size_tup" = size_tup, "size_uct" = x, "sdeps" = sdeps, "n0" = NULL, "nuct" = NULL))
        varsoln_size_uct <- var_bench(soln_size_uct)
        list("x" = x, "soln" = soln_size_uct, "varsoln" = varsoln_size_uct)
      })
    }
    if(plotyvar == "n") {
      df <- lapply(v, function(vi) {
        data.frame(x = vi$x, group = c("n0", "n1", "n2"),
                   y = c(vi[["soln"]]$n0, vi[["soln"]]$nuct, 1 - (vi[["soln"]]$n0 +  vi[["soln"]]$nuct)))
      }) %>% bind_rows()
    } else if(plotyvar == "var") {
      df <- lapply(v, function(vi) {
        data.frame(x = vi$x, group = c("var0", "var1", "var2"),
                   y = c(vi[["varsoln"]]["var_betatup"], vi[["varsoln"]]["var_betatupce"], vi[["varsoln"]]["var_betauct"]))
      }) %>% bind_rows()
    }
    if(plotyvar == "n") {
      labels <- c(expression(N[Control]/N), expression(N[UCT]/N), expression(N[Program]/N))
    } else if(plotyvar == "var") {
      # labels <- c(expression(paste("Var[", hat("\u03B2")[2], "]")),
      #             expression(paste("Var[", hat("\u03B2")[2], " \u2212 ", over("Size"["Program"], "Size"["UCT"]), "", hat("\u03B2")[1], "]")))
      labels <- c(expression(paste("Var[", hat("\u03B2'")[2], "]")),
                  expression(paste("Var[", hat("\u03B4'")[2], "]")),
                  expression(paste("Var[", hat("\u03B2'")[1], "]")))
    }
    #### FIX BELOW HERE ####
    # GENERATE PLOT
    g <- ggplot(df, aes(x = x, y = y, group = group, col = group)) +
      geom_line() +
      geom_vline(xintercept = case_when(plotxvar == "omega1" ~ omega1, plotxvar == "omega2" ~ omega2,
                                        plotxvar == "size_uct" ~ size_uct), linetype = 2) +
      xlim(0, case_when(plotxvar == "omega1" ~ 1, plotxvar == "omega2" ~ 1,
                        plotxvar == "size_uct" ~ 2*max(size_uct, size_uct_other, 2))) +
      # ylim(0, case_when(plotyvar == "n" ~ 1, plotyvar == "var" ~ 2 * max(varsoln, varsoln_other))) +
      xlab(case_when(plotxvar == "omega1" ~ expression("\u03C9"[1]), plotxvar == "omega2" ~ expression("\u03C9"[2]),
                     plotxvar == "size_uct" ~ expression(Size[UCT]))) +
      # ylab(case_when(plotyvar == "n" ~ "Fraction of obs.", plotyvar == "var" ~ "Variance")) +
      scale_color_discrete(labels = labels) +
      theme_bw(base_size = 16) +
      theme(legend.title = element_blank(), legend.position = c(0.7, 0.8)) +
      ggtitle(paste0("Panel ", j))
    if(plotyvar == "n") {
      g <- g + ylim(0, 1) + ylab("Fraction of obs.")
    } else if(plotyvar == "var") {
      g <- g + scale_y_log10(limits = c(0.8 * min(df$y), 2 * max(varsoln, varsoln_other)), name = "Variance (log scale)")
    }
    return(g)
  }
  
  mz_genpaneltable <- function(isize_uct1, isize_uct2, iomega1, iomega2) {
    param1 <- as.numeric(input$mz_param1)
    size_uct1 <- as.numeric(isize_uct1)
    size_uct2 <- as.numeric(isize_uct2)
    omega1 <- as.numeric(iomega1)
    omega2 <- as.numeric(iomega2)
    size_tupmin <- as.numeric(input$mz_size_tupmin)
    size_tupmax <- as.numeric(input$mz_size_tupmax)
    sdeps <- as.numeric(input$mz_sdeps)
    soln <- opt_mz18(list("size_tupmin" = size_tupmin, "size_tupmax" = size_tupmax, "omega1" = omega1, "omega2" = omega2, "rho" = param1, "size_uct1" = size_uct1,
                          "size_uct2" = size_uct2, "sdeps" = sdeps, "n0" = NULL, "nuct1" = NULL, "nuct2" = NULL))
    varsoln <- var_mz18(soln)
    tab <- data.frame(c("Fraction Control" = soln$n0, "Fraction Small UCT" = soln$nuct1, "Fraction Large UCT" = soln$nuct2,
                        "Fraction Program" = 1 - (soln$n0 + soln$nuct1 + soln$nuct2),
                        "Var[Est. Program Effect]" = varsoln["var_betatup"] %>% unname,
                        "Var[Est. Program Effect Rel. to Cash Benchmark]" = varsoln["var_betatupce"] %>% unname,
                        "Var[Large UCT Marg. Effect Rel. to Small UCT]" = varsoln["var_betareluct"] %>% unname)) %>%
      as.matrix
    colnames(tab) <- " "
    # paste0("\\[", print.xtableFtable(xtableFtable(tab, method = "compact"), booktabs = T), "\\]")
    tab
  }
  
  mz_genpowertable <- function(powervar, isize_uct1, isize_uct2, iomega1, iomega2) {
    param1 <- as.numeric(input$mz_param1)
    size_uct1 <- as.numeric(isize_uct1)
    size_uct2 <- as.numeric(isize_uct2)
    omega1 <- as.numeric(iomega1)
    omega2 <- as.numeric(iomega2)
    size_tupmin <- as.numeric(input$mz_size_tupmin)
    size_tupmax <- as.numeric(input$mz_size_tupmax)
    sdeps <- as.numeric(input$mz_sdeps)
    soln <- opt_mz18(list("size_tupmin" = size_tupmin, "size_tupmax" = size_tupmax, "omega1" = omega1, "omega2" = omega2, "rho" = param1, "size_uct1" = size_uct1,
                          "size_uct2" = size_uct2, "sdeps" = sdeps, "n0" = NULL, "nuct1" = NULL, "nuct2" = NULL))
    varsoln <- var_mz18(soln)
    alpha <- as.numeric(input$mz_alpha)
    power <- as.numeric(input$mz_power)
    N <- as.numeric(input$mz_samplesize)
    mde <- c(as.numeric(input$mz_mde_betatup), as.numeric(input$mz_mde_betatupce), as.numeric(input$mz_mde_betareluct))
    if(powervar == "mde") {
      out <- unname((qnorm(1 - alpha/2) + qnorm(power)) * sqrt(varsoln / N))
      # outnames <- c("MDE (Program Effect)", "MDE (Program Effect Rel. to Cash Benchmark)", "MDE (Large UCT Marg. Effect Rel. to Small UCT)")
    } else if(powervar == "n") {
      out <- unname((qnorm(1 - alpha/2) + qnorm(power))^2 * varsoln / (mde^2))
    }
    outnames <- c("Program Effect", "Program Effect Rel. to Cash Benchmark", "Large UCT Marg. Effect Rel. to Small UCT")
    tab <- data.frame(setNames(out, outnames)) %>% as.matrix
    colnames(tab) <- " "
    tab
  }
  
  mz_genplot <- function(j, plotxvar, plotyvar, isize_uct1, isize_uct2, iomega1, iomega2, isize_uct1_other, isize_uct2_other, iomega1_other, iomega2_other) {
    param1 <- as.numeric(input$mz_param1)
    size_tupmin <- as.numeric(input$mz_size_tupmin)
    size_tupmax <- as.numeric(input$mz_size_tupmax)
    sdeps <- as.numeric(input$mz_sdeps)
    size_uct1 <- as.numeric(isize_uct1)
    size_uct2 <- as.numeric(isize_uct2)
    omega1 <- as.numeric(iomega1)
    omega2 <- as.numeric(iomega2)
    size_uct1_other <- as.numeric(isize_uct1_other)
    size_uct2_other <- as.numeric(isize_uct2_other)
    omega1_other <- as.numeric(iomega1_other)
    omega2_other <- as.numeric(iomega2_other)
    soln <- opt_mz18(list("size_tupmin" = size_tupmin, "size_tupmax" = size_tupmax, "omega1" = omega1, "omega2" = omega2, "rho" = param1, "size_uct1" = size_uct1,
                          "size_uct2" = size_uct2, "sdeps" = sdeps, "n0" = NULL, "nuct1" = NULL, "nuct2" = NULL))
    varsoln <- var_mz18(soln)
    soln_other <- opt_mz18(list("size_tupmin" = size_tupmin, "size_tupmax" = size_tupmax, "omega1" = omega1_other, "omega2" = omega2_other, "rho" = param1, "size_uct1" = size_uct1_other,
                                "size_uct2" = size_uct2_other, "sdeps" = sdeps, "n0" = NULL, "nuct1" = NULL, "nuct2" = NULL))
    varsoln_other <- var_mz18(soln_other)
    # GENERATE PLOT DATA
    if(plotxvar == "omega1") {
      v <- lapply(seq(0.0125, 0.9875, 0.025), function(x) {
        soln_omega1 <- opt_mz18(list("size_tupmin" = size_tupmin, "size_tupmax" = size_tupmax, "omega1" = x, "omega2" = (1 - x) * omega2 / (1 - omega1),
                                     "rho" = param1, "size_uct1" = size_uct1,
                                     "size_uct2" = size_uct2, "sdeps" = sdeps, "n0" = NULL, "nuct1" = NULL, "nuct2" = NULL))
        varsoln_omega1 <- var_mz18(soln_omega1)
        list("x" = x, "soln" = soln_omega1, "varsoln" = varsoln_omega1)
      })
    } else if(plotxvar == "omega2") {
      v <- lapply(seq(0.0125, 0.9875, 0.025), function(x) {
        soln_omega2 <- opt_mz18(list("size_tupmin" = size_tupmin, "size_tupmax" = size_tupmax, "omega1" = (1 - x) * omega1 / (1 - omega2), "omega2" = x,
                                     "rho" = param1, "size_uct1" = size_uct1,
                                     "size_uct2" = size_uct2, "sdeps" = sdeps, "n0" = NULL, "nuct1" = NULL, "nuct2" = NULL))
        varsoln_omega2 <- var_mz18(soln_omega2)
        list("x" = x, "soln" = soln_omega2, "varsoln" = varsoln_omega2)
      })
    } else if(plotxvar == "size_uct1") {
      v <- lapply(seq(0, size_uct2, length.out = 41), function(x) {
        soln_size_uct1 <- opt_mz18(list("size_tupmin" = size_tupmin, "size_tupmax" = size_tupmax, "omega1" = omega1, "omega2" = omega2, "rho" = param1, "size_uct1" = x,
                                        "size_uct2" = size_uct2, "sdeps" = sdeps, "n0" = NULL, "nuct1" = NULL, "nuct2" = NULL))
        varsoln_size_uct1 <- var_mz18(soln_size_uct1)
        list("x" = x, "soln" = soln_size_uct1, "varsoln" = varsoln_size_uct1)
      })
    } else if(plotxvar == "size_uct2") {
      v <- lapply(seq(size_uct1, 2*max(2, size_uct2, size_uct2_other), length.out = 41), function(x) {
        soln_size_uct2 <- opt_mz18(list("size_tupmin" = size_tupmin, "size_tupmax" = size_tupmax, "omega1" = omega1, "omega2" = omega2, "rho" = param1, "size_uct1" = size_uct1,
                                        "size_uct2" = x, "sdeps" = sdeps, "n0" = NULL, "nuct1" = NULL, "nuct2" = NULL))
        varsoln_size_uct2 <- var_mz18(soln_size_uct2)
        list("x" = x, "soln" = soln_size_uct2, "varsoln" = varsoln_size_uct2)
      })
    }
    if(plotyvar == "n") {
      df <- lapply(v, function(vi) {
        data.frame(x = vi$x, group = c("n0", "n1", "n2", "n3"),
                   y = c(vi[["soln"]]$n0, vi[["soln"]]$nuct1, vi[["soln"]]$nuct2, 1 - (vi[["soln"]]$n0 +  vi[["soln"]]$nuct1 + vi[["soln"]]$nuct2)))
      }) %>% bind_rows()
    } else if(plotyvar == "var") {
      df <- lapply(v, function(vi) {
        data.frame(x = vi$x, group = c("var0", "var1", "var2"),
                   y = c(vi[["varsoln"]]["var_betatup"], vi[["varsoln"]]["var_betatupce"], vi[["varsoln"]]["var_betareluct"]))
      }) %>% bind_rows()
    }
    if(plotyvar == "n") {
      labels <- c(expression(N[Control]/N), expression(N[Small~UCT]/N), expression(N[Large~UCT]/N), expression(N[Program]/N))
    } else if(plotyvar == "var") {
      labels <- c(expression(paste("Var[", hat("\u03B2'")[3], "]")),
                  expression(paste("Var[", hat("\u03B4'")[3], "]")),
                  #expression(paste("Var[", over(hat("\u03B2")[2], "Size"["Large UCT"]),  " \u2212 ", over(hat("\u03B2")[1], "Size"["Small UCT"]), "]")))
                  expression(paste("Var[", hat("\u03B2'")[2],  " \u2212 ", hat("\u03B2'")[1], "]")))
    }
    # GENERATE PLOT
    g <- ggplot(df, aes(x = x, y = y, group = group, col = group)) +
      geom_line() +
      geom_vline(xintercept = case_when(plotxvar == "omega1" ~ omega1, plotxvar == "omega2" ~ omega2,
                                        plotxvar == "size_uct1" ~ size_uct1, plotxvar == "size_uct2" ~ size_uct2), linetype = 2) +
      xlim(0, case_when(plotxvar == "omega1" ~ 1, plotxvar == "omega2" ~ 1, plotxvar == "size_uct1" ~ 2*max(size_uct2, size_uct2_other, 2),
                        plotxvar == "size_uct2" ~ 2*max(size_uct2, size_uct2_other, 2))) +
      # ylim(0, case_when(plotyvar == "n" ~ 1, plotyvar == "var" ~ 2 * max(varsoln, varsoln_other))) +
      xlab(case_when(plotxvar == "omega1" ~ expression("\u03C9"[1]), plotxvar == "omega2" ~ expression("\u03C9"[2]),
                     plotxvar == "size_uct1" ~ expression(Size[Small~UCT]), plotxvar == "size_uct2" ~ expression(Size[Large~UCT]))) +
      # ylab(case_when(plotyvar == "n" ~ "Fraction of obs.", plotyvar == "var" ~ "Variance")) +
      scale_color_discrete(labels = labels) +
      theme_bw(base_size = 16) +
      theme(legend.title = element_blank(), legend.position = c(0.7, 0.8)) +
      ggtitle(paste0("Panel ", j))
    if(plotyvar == "n") {
      g <- g + ylim(0, 1) + ylab("Fraction of obs.")
    } else if(plotyvar == "var") {
      g <- g + scale_y_log10(limits = c(0.8 * min(df$y), 2 * max(varsoln, varsoln_other)), name = "Variance (log scale)")
    }
    return(g)
  }
  
  output$nl_panel1table <- renderTable({ nl_genpaneltable(input$nl_panel1_size_uct1, input$nl_panel1_size_uct2) }, rownames = T, digits = 4)
  output$nl_panel2table <- renderTable({ nl_genpaneltable(input$nl_panel2_size_uct1, input$nl_panel2_size_uct2) }, rownames = T, digits = 4)
  output$bm_panel1table <- renderTable({ bm_genpaneltable(input$bm_panel1_size_uct, input$bm_panel1_omega1, input$bm_panel1_omega2) }, rownames = T, digits = 4)
  output$bm_panel2table <- renderTable({ bm_genpaneltable(input$bm_panel2_size_uct, input$bm_panel2_omega1, input$bm_panel2_omega2) }, rownames = T, digits = 4)
  output$mz_panel1table <- renderTable({ mz_genpaneltable(input$mz_panel1_size_uct1, input$mz_panel1_size_uct2,
                                                          input$mz_panel1_omega1, input$mz_panel1_omega2) }, rownames = T, digits = 4)
  output$mz_panel2table <- renderTable({ mz_genpaneltable(input$mz_panel2_size_uct1, input$mz_panel2_size_uct2,
                                                          input$mz_panel2_omega1, input$mz_panel2_omega2) }, rownames = T, digits = 4)
  
  output$nl_panel1_table_mde <- renderTable({ nl_genpowertable("mde", input$nl_panel1_size_uct1, input$nl_panel1_size_uct2) }, rownames = T, digits = 4)
  output$nl_panel1_table_samplesize <- renderTable({ nl_genpowertable("n", input$nl_panel1_size_uct1, input$nl_panel1_size_uct2) }, rownames = T, digits = 0)
  output$bm_panel1_table_mde <- renderTable({ bm_genpowertable("mde", input$bm_panel1_size_uct, input$bm_panel1_omega1, input$bm_panel1_omega2) }, rownames = T, digits = 4)
  output$bm_panel1_table_samplesize <- renderTable({ bm_genpowertable("n", input$bm_panel1_size_uct, input$bm_panel1_omega1, input$bm_panel1_omega2) }, rownames = T, digits = 0)
  output$mz_panel1_table_mde <- renderTable({ mz_genpowertable("mde", input$mz_panel1_size_uct1, input$mz_panel1_size_uct2,
                                                               input$mz_panel1_omega1, input$mz_panel1_omega2) }, rownames = T, digits = 4)
  output$mz_panel1_table_samplesize <- renderTable({ mz_genpowertable("n", input$mz_panel1_size_uct1, input$mz_panel1_size_uct2,
                                                                      input$mz_panel1_omega1, input$mz_panel1_omega2) }, rownames = T, digits = 0)
  
  output$nl_panel1_plot_n_size_uct1 <- renderPlot({ nl_genplot(1, "size_uct1", "n", input$nl_panel1_size_uct1, input$nl_panel1_size_uct2,
                                                               input$nl_panel2_size_uct1, input$nl_panel2_size_uct2) })
  output$nl_panel1_plot_n_size_uct2 <- renderPlot({ nl_genplot(1, "size_uct2", "n", input$nl_panel1_size_uct1, input$nl_panel1_size_uct2,
                                                               input$nl_panel2_size_uct1, input$nl_panel2_size_uct2) })
  output$nl_panel1_plot_var_size_uct1 <- renderPlot({ nl_genplot(1, "size_uct1", "var", input$nl_panel1_size_uct1, input$nl_panel1_size_uct2,
                                                                 input$nl_panel2_size_uct1, input$nl_panel2_size_uct2) })
  output$nl_panel1_plot_var_size_uct2 <- renderPlot({ nl_genplot(1, "size_uct2", "var", input$nl_panel1_size_uct1, input$nl_panel1_size_uct2,
                                                                 input$nl_panel2_size_uct1, input$nl_panel2_size_uct2) })
  output$nl_panel2_plot_n_size_uct1 <- renderPlot({ nl_genplot(2, "size_uct1", "n", input$nl_panel2_size_uct1, input$nl_panel2_size_uct2,
                                                               input$nl_panel1_size_uct1, input$nl_panel1_size_uct2) })
  output$nl_panel2_plot_n_size_uct2 <- renderPlot({ nl_genplot(2, "size_uct2", "n", input$nl_panel2_size_uct1, input$nl_panel2_size_uct2,
                                                               input$nl_panel1_size_uct1, input$nl_panel1_size_uct2) })
  output$nl_panel2_plot_var_size_uct1 <- renderPlot({ nl_genplot(2, "size_uct1", "var", input$nl_panel2_size_uct1, input$nl_panel2_size_uct2,
                                                                 input$nl_panel1_size_uct1, input$nl_panel1_size_uct2) })
  output$nl_panel2_plot_var_size_uct2 <- renderPlot({ nl_genplot(2, "size_uct2", "var", input$nl_panel2_size_uct1, input$nl_panel2_size_uct2,
                                                                 input$nl_panel1_size_uct1, input$nl_panel1_size_uct2) })
  
  output$bm_panel1_plot_n_size_uct <- renderPlot({ bm_genplot(1, "size_uct", "n", input$bm_panel1_size_uct, input$bm_panel1_omega1, input$bm_panel1_omega2,
                                                              input$bm_panel2_size_uct, input$bm_panel2_omega1, input$bm_panel2_omega2) })
  output$bm_panel1_plot_n_omega1 <- renderPlot({ bm_genplot(1, "omega1", "n", input$bm_panel1_size_uct, input$bm_panel1_omega1, input$bm_panel1_omega2,
                                                           input$bm_panel2_size_uct, input$bm_panel2_omega1, input$bm_panel2_omega2) })
  output$bm_panel1_plot_n_omega2 <- renderPlot({ bm_genplot(1, "omega2", "n", input$bm_panel1_size_uct, input$bm_panel1_omega1, input$bm_panel1_omega2,
                                                           input$bm_panel2_size_uct, input$bm_panel2_omega1, input$bm_panel2_omega2) })
  output$bm_panel1_plot_var_size_uct <- renderPlot({ bm_genplot(1, "size_uct", "var", input$bm_panel1_size_uct, input$bm_panel1_omega1, input$bm_panel1_omega2,
                                                                input$bm_panel2_size_uct, input$bm_panel2_omega1, input$bm_panel2_omega2) })
  output$bm_panel1_plot_var_omega1 <- renderPlot({ bm_genplot(1, "omega1", "var", input$bm_panel1_size_uct, input$bm_panel1_omega1, input$bm_panel1_omega2,
                                                             input$bm_panel2_size_uct, input$bm_panel2_omega1, input$bm_panel2_omega2) })
  output$bm_panel1_plot_var_omega2 <- renderPlot({ bm_genplot(1, "omega2", "var", input$bm_panel1_size_uct, input$bm_panel1_omega1, input$bm_panel1_omega2,
                                                             input$bm_panel2_size_uct, input$bm_panel2_omega1, input$bm_panel2_omega2) })
  output$bm_panel2_plot_n_size_uct <- renderPlot({ bm_genplot(2, "size_uct", "n", input$bm_panel2_size_uct, input$bm_panel2_omega1, input$bm_panel2_omega2,
                                                              input$bm_panel1_size_uct, input$bm_panel1_omega1, input$bm_panel1_omega2) })
  output$bm_panel2_plot_n_omega1 <- renderPlot({ bm_genplot(2, "omega1", "n", input$bm_panel2_size_uct, input$bm_panel2_omega1, input$bm_panel2_omega2,
                                                           input$bm_panel1_size_uct, input$bm_panel1_omega1, input$bm_panel1_omega2) })
  output$bm_panel2_plot_n_omega2 <- renderPlot({ bm_genplot(2, "omega2", "n", input$bm_panel2_size_uct, input$bm_panel2_omega1, input$bm_panel2_omega2,
                                                           input$bm_panel1_size_uct, input$bm_panel1_omega1, input$bm_panel1_omega2) })
  output$bm_panel2_plot_var_size_uct <- renderPlot({ bm_genplot(2, "size_uct", "var", input$bm_panel2_size_uct, input$bm_panel2_omega1, input$bm_panel2_omega2,
                                                                input$bm_panel1_size_uct, input$bm_panel1_omega1, input$bm_panel1_omega2) })
  output$bm_panel2_plot_var_omega1 <- renderPlot({ bm_genplot(2, "omega1", "var", input$bm_panel2_size_uct, input$bm_panel2_omega1, input$bm_panel2_omega2,
                                                             input$bm_panel1_size_uct, input$bm_panel1_omega1, input$bm_panel1_omega2) })
  output$bm_panel2_plot_var_omega2 <- renderPlot({ bm_genplot(2, "omega2", "var", input$bm_panel2_size_uct, input$bm_panel2_omega1, input$bm_panel2_omega2,
                                                             input$bm_panel1_size_uct, input$bm_panel1_omega1, input$bm_panel1_omega2) })
  
  output$mz_panel1_plot_n_size_uct1 <- renderPlot({ mz_genplot(1, "size_uct1", "n",
                                                               input$mz_panel1_size_uct1, input$mz_panel1_size_uct2, input$mz_panel1_omega1, input$mz_panel1_omega2,
                                                               input$mz_panel2_size_uct1, input$mz_panel2_size_uct2, input$mz_panel2_omega1, input$mz_panel2_omega2) })
  output$mz_panel1_plot_n_size_uct2 <- renderPlot({ mz_genplot(1, "size_uct2", "n",
                                                               input$mz_panel1_size_uct1, input$mz_panel1_size_uct2, input$mz_panel1_omega1, input$mz_panel1_omega2,
                                                               input$mz_panel2_size_uct1, input$mz_panel2_size_uct2, input$mz_panel2_omega1, input$mz_panel2_omega2) })
  output$mz_panel1_plot_n_omega1 <- renderPlot({ mz_genplot(1, "omega1", "n",
                                                            input$mz_panel1_size_uct1, input$mz_panel1_size_uct2, input$mz_panel1_omega1, input$mz_panel1_omega2,
                                                            input$mz_panel2_size_uct1, input$mz_panel2_size_uct2, input$mz_panel2_omega1, input$mz_panel2_omega2) })
  output$mz_panel1_plot_n_omega2 <- renderPlot({ mz_genplot(1, "omega2", "n",
                                                            input$mz_panel1_size_uct1, input$mz_panel1_size_uct2, input$mz_panel1_omega1, input$mz_panel1_omega2,
                                                            input$mz_panel2_size_uct1, input$mz_panel2_size_uct2, input$mz_panel2_omega1, input$mz_panel2_omega2) })
  output$mz_panel1_plot_var_size_uct1 <- renderPlot({ mz_genplot(1, "size_uct1", "var",
                                                                 input$mz_panel1_size_uct1, input$mz_panel1_size_uct2, input$mz_panel1_omega1, input$mz_panel1_omega2,
                                                                 input$mz_panel2_size_uct1, input$mz_panel2_size_uct2, input$mz_panel2_omega1, input$mz_panel2_omega2) })
  output$mz_panel1_plot_var_size_uct2 <- renderPlot({ mz_genplot(1, "size_uct2", "var",
                                                                 input$mz_panel1_size_uct1, input$mz_panel1_size_uct2, input$mz_panel1_omega1, input$mz_panel1_omega2,
                                                                 input$mz_panel2_size_uct1, input$mz_panel2_size_uct2, input$mz_panel2_omega1, input$mz_panel2_omega2) })
  output$mz_panel1_plot_var_omega1 <- renderPlot({ mz_genplot(1, "omega1", "var",
                                                              input$mz_panel1_size_uct1, input$mz_panel1_size_uct2, input$mz_panel1_omega1, input$mz_panel1_omega2,
                                                              input$mz_panel2_size_uct1, input$mz_panel2_size_uct2, input$mz_panel2_omega1, input$mz_panel2_omega2) })
  output$mz_panel1_plot_var_omega2 <- renderPlot({ mz_genplot(1, "omega2", "var",
                                                              input$mz_panel1_size_uct1, input$mz_panel1_size_uct2, input$mz_panel1_omega1, input$mz_panel1_omega2,
                                                              input$mz_panel2_size_uct1, input$mz_panel2_size_uct2, input$mz_panel2_omega1, input$mz_panel2_omega2) })
  output$mz_panel2_plot_n_size_uct1 <- renderPlot({ mz_genplot(2, "size_uct1", "n",
                                                               input$mz_panel2_size_uct1, input$mz_panel2_size_uct2, input$mz_panel2_omega1, input$mz_panel2_omega2,
                                                               input$mz_panel1_size_uct1, input$mz_panel1_size_uct2, input$mz_panel1_omega1, input$mz_panel1_omega2) })
  output$mz_panel2_plot_n_size_uct2 <- renderPlot({ mz_genplot(2, "size_uct2", "n",
                                                               input$mz_panel2_size_uct1, input$mz_panel2_size_uct2, input$mz_panel2_omega1, input$mz_panel2_omega2,
                                                               input$mz_panel1_size_uct1, input$mz_panel1_size_uct2, input$mz_panel1_omega1, input$mz_panel1_omega2) })
  output$mz_panel2_plot_n_omega1 <- renderPlot({ mz_genplot(2, "omega1", "n",
                                                            input$mz_panel2_size_uct1, input$mz_panel2_size_uct2, input$mz_panel2_omega1, input$mz_panel2_omega2,
                                                            input$mz_panel1_size_uct1, input$mz_panel1_size_uct2, input$mz_panel1_omega1, input$mz_panel1_omega2) })
  output$mz_panel2_plot_n_omega2 <- renderPlot({ mz_genplot(2, "omega2", "n",
                                                            input$mz_panel2_size_uct1, input$mz_panel2_size_uct2, input$mz_panel2_omega1, input$mz_panel2_omega2,
                                                            input$mz_panel1_size_uct1, input$mz_panel1_size_uct2, input$mz_panel1_omega1, input$mz_panel1_omega2) })
  output$mz_panel2_plot_var_size_uct1 <- renderPlot({ mz_genplot(2, "size_uct1", "var",
                                                                 input$mz_panel2_size_uct1, input$mz_panel2_size_uct2, input$mz_panel2_omega1, input$mz_panel2_omega2,
                                                                 input$mz_panel1_size_uct1, input$mz_panel1_size_uct2, input$mz_panel1_omega1, input$mz_panel1_omega2) })
  output$mz_panel2_plot_var_size_uct2 <- renderPlot({ mz_genplot(2, "size_uct2", "var",
                                                                 input$mz_panel2_size_uct1, input$mz_panel2_size_uct2, input$mz_panel2_omega1, input$mz_panel2_omega2,
                                                                 input$mz_panel1_size_uct1, input$mz_panel1_size_uct2, input$mz_panel1_omega1, input$mz_panel1_omega2) })
  output$mz_panel2_plot_var_omega1 <- renderPlot({ mz_genplot(2, "omega1", "var",
                                                              input$mz_panel2_size_uct1, input$mz_panel2_size_uct2, input$mz_panel2_omega1, input$mz_panel2_omega2,
                                                              input$mz_panel1_size_uct1, input$mz_panel1_size_uct2, input$mz_panel1_omega1, input$mz_panel1_omega2) })
  output$mz_panel2_plot_var_omega2 <- renderPlot({ mz_genplot(2, "omega2", "var",
                                                              input$mz_panel2_size_uct1, input$mz_panel2_size_uct2, input$mz_panel2_omega1, input$mz_panel2_omega2,
                                                              input$mz_panel1_size_uct1, input$mz_panel1_size_uct2, input$mz_panel1_omega1, input$mz_panel1_omega2) })
}

shinyApp(ui, server)
