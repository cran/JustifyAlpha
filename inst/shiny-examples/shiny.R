library(shiny)
library(pwr)
library(ggplot2)
library(shinydashboard)
library(BayesFactor)
library(JustifyAlpha)

bf_bic <- function(Fval, df1, df2, repeated=FALSE, report.as="BF10") {
  if (repeated==FALSE) {
    N = df1+df2+1
  }
  else {
    N = df1+df2
  }
  
  bf = sqrt(N^df1*(1+Fval*df1/df2)^(-1*N))
  
  if (report.as=="BF01"){
    return(c(B01=bf))
  }
  else {
    return(c(B10=1/bf))
  }
}

ui <- dashboardPage(
  dashboardHeader(title = "Justify Your Alpha"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Balance/Minimize Error Rates", tabName = "error",
               menuSubItem("Justify Alpha", tabName = "minimize"),
               menuSubItem("Justify Sample Size and Alpha", tabName = "sample"),
               icon = icon("calculator")),
      menuItem("N as a Function of Sample Size", tabName = "alpha_sample_size",
               menuSubItem("t-test", tabName = "ttest"),
               menuSubItem("ANOVA", tabName = "anova"),
               icon = icon("calculator")),
      menuItem("About", tabName = "about", icon = icon("info"))
    )
  ),
  dashboardBody(
    tabItems(
      # First tab content
      tabItem(tabName = "minimize",
              fluidRow(
                box(
                  title = "Input parameters and press 'Calculate'",
                  selectInput("error", "Minimize or Balance Error Rates?:",
                              c("Balance" = "balance",
                                "Minimize" = "minimize"
                              )),
                  numericInput("costT1T2", "Relative cost Type 1 and Type 2 errors:", 4, min = 0),
                  textOutput("outT1T2"),
                  br(),
                  numericInput("priorH1H0", "Prior Odds of H1 compared to H0:", 1, min = 0),
                  textOutput("outH1H0"),
                  br(),
                  textAreaInput("power_function", "Power function:", "pwr::pwr.t.test(d = 0.5, n = 64, sig.level = x, type = 'two.sample', alternative = 'two.sided')$power", width = '400px', height = '200px'),
                  actionButton("power_start", "Calculate")
                ),
                infoBoxOutput("alpha1Box"),
                infoBoxOutput("beta1Box"),
                infoBoxOutput("error1Box"),
                box(plotOutput("plot1")),
                box(title = "Explanation",
                    status = "warning",
                    solidHeader = TRUE, collapsible = TRUE,
                    "This Shiny app accompanies Maier & Lakens (2021). Justify Your Alpha: A Primer on Two Practical Approaches. For a full explanation on how to use this software, read the paper or the vignettes."),
                box(title = "Power functions",
                    status = "info",
                    solidHeader = TRUE, collapsible = TRUE,
                    "The trickiest thing of using this Shiny app is entering the correct power function. You can provide an analytic power function, either programmed yourself, or from an existing package loading on the server. Then, make sure the alpha value is not set, but specified as x, and that the function itself returns a single value, the power of the test. Finally, if you use existing power functions the shiny app needs to know which package this function is from, and thus the call to the function needs to be precended by the package and '::', so 'pwr::' or 'TOSTER::'. Some examples that work are provided below.", tags$br(), tags$br(),
                    "pwr::pwr.anova.test(n = 100, k = 2, f = 0.171875, sig.level = x)$power", tags$br(), tags$br(),
                    "TOSTER::powerTOSTtwo(alpha=x, N=200, low_eqbound_d=-0.4, high_eqbound_d=0.4)", tags$br(), tags$br(),
                    "For a more challenging power function, we can use the Superpower package by Daniel Lakens and Aaron Caldwell. The power function in the ANOVAexact function is based on a simulation, which takes a while to perform. The optimization function used in this Shiny app needs to perform the power calculation multiple times. Thus, the result takes a minutes to calculate. Press calculate, and check the results 5 to 10 minutes later. Furthermore, the output of the ANOVA_exact function prints power as 80%, not 0.8, and thus we actually have to divide the power value by 100 for the Shiny app to return the correct results. Nevertheless, it works if you are very patient.", tags$br(), tags$br(),
                    "Superpower::ANOVA_exact( (Superpower::ANOVA_design(design = '2b', n = 64, mu = c(0, 0.5), sd = 1, plot = FALSE)), alpha_level = x, verbose = FALSE)$main_results$power/100"
                )
              )
      ),
      # Second tab content
      tabItem(tabName = "sample",
              fluidRow(
                box(
                  title = "Input parameters and press 'Calculate'",
                  selectInput("error2", "Minimize or Balance Error Rates?:",
                              c("Balance" = "balance",
                                "Minimize" = "minimize"
                              )),
                  numericInput("costT1T22", "Relative cost Type 1 and Type 2 errors:", 4, min = 0),
                  textOutput("outT1T22"),
                  br(),
                  numericInput("priorH1H02", "Prior Odds of H1 compared to H0:", 1, min = 0),
                  textOutput("outH1H02"),
                  br(),
                  numericInput("errorrate2", "Desired Weighted Combined Error Rate", 0.05, min = 0, max = 1),
                  textAreaInput("power_function2", "Power function:", "pwr::pwr.t.test(d = 0.5, n = sample_n, sig.level = x, type = 'two.sample', alternative = 'two.sided')$power", width = '400px', height = '200px'),
                  actionButton("power_start2", "Calculate")
                ),
                infoBoxOutput("alpha1Box2"),
                infoBoxOutput("beta1Box2"),
                infoBoxOutput("sampleBox2"),
                infoBoxOutput("error1Box2"),
                box(plotOutput("plot2")),
                box(title = "Explanation",
                    status = "warning",
                    solidHeader = TRUE, collapsible = TRUE,
                    "Cohen (1988) considered a Type 1 error rate of 5% and a Type 2 error rate of 20% balanced. The reason for this was that instead of weighing both types of errors equally, he felt 'Type I errors are of the order of four times as serious as Type II errors.' This situation is illustrated in the default settings of the app. If the cost of a Type 1 error is 4 times as large as the cost of a Type 2 error, and we collect 64 participants in each condition of a two-sided t-test, that alpha is 0.05 and power is 0.80.", tags$br(), tags$br(),
                    "If we design 2000 studies like this, the number of Type 1 and Type 2 errors we make depend on how often the null hypothesis is true, and how often the alternative hypothesis is true. Let's assume both are equally likely for now. This means that in 1000 studies the null hypothesis is true, and we will make 50 Type 1 errors. In 1000 studies the alternative hypothesis is true, and we will make 100%-80% = 20% Type 2 errors, so in 200 studies we will not observe a significant result even if there is a true effect. Combining Type 1 and Type 2 errors, in the long run, we should expect 250 of our 2000 studies to yield an error."
                ),
                box(title = "Power functions",
                    status = "info",
                    solidHeader = TRUE, collapsible = TRUE,
                    "The trickiest thing of using this Shiny app is entering the correct power function. You can provide an analytic power function, either programmed yourself, or from an existing package loading on the server. Then, make sure the alpha value is not set, but specified as x, and that the sample is not set but specified as 'sample_n'. In additino, make sure that the function itself returns a single value, the power of the test. Finally, if you use existing power functions the shiny app needs to know which package this function is from, and thus the call to the function needs to be precended by the package and '::', so 'pwr::' or 'TOSTER::' or 'ANOVApower::'. Some examples that work are provided below.", tags$br(), tags$br(),
                    "TOSTER::powerTOSTtwo(alpha=x, N=200, low_eqbound_d=-0.4, high_eqbound_d=0.4)", tags$br(), tags$br(),
                    "pwr::pwr.anova.test(n = sample_n, k = 2, f = 0.171875, sig.level = x)$power", tags$br(), tags$br(),
                )
              )
      ),
     tabItem(tabName = "ttest",
              fluidRow(
                box(
                  title = "Input parameters",
                  numericInput("n1", "n1", 50, min = 0),
                  numericInput("n2", "n2 (0 for one sample/paired test)", 0, min = 0),
                  selectInput("one.sided", "One sided or two sided?",  c("two sided" = FALSE,
                                                                    "one sided" = TRUE)),
                  sliderInput("evidence", "How much more likely should the data at least be under the alternative hypothesis?",
                              min = 1,
                              max = 10,
                              value = 3),
                  numericInput("rscale", "Cauchy scale (advanced)", 0.707, min = 0.1),
                  textOutput("likelyttest"),
                  br(),
                  actionButton("power_start3", "Calculate")
                  ),
                infoBoxOutput("ttestbox"),
                box(plotOutput("plotttest")),
                #box(plotOutput("ttestplot")),
                box(title = "Explanation",
                    status = "warning",
                    solidHeader = TRUE, collapsible = TRUE,
                    "The idea behind this recommendation is discussed most extensively by Leamer, 1978. He writes 'The rule of thumb quite popular now, that is, setting the significance level arbitrarily to .05, is shown to be deficient in the sense that from every reasonable viewpoint the significance level should be a decreasing function of sample size.' This was already recognized by Jeffreys (1939), who discusses ways to set the alpha level in Neyman-Pearson statistics: 'We should therefore get the best result, with any distribution of alpha, by some form that makes the ratio of the critical value to the standard error increase with n. It appears then that whatever the distribution may be, the use of a fixed P limit cannot be the one that will make the smallest number of mistakes.'", tags$br(), tags$br(),
                    "The goal is to prevent Lindley's paradox (https://en.wikipedia.org/wiki/Lindley%27s_paradox). This is explained in more detail in week 1 of Daniel's MOOC (https://www.coursera.org/learn/statistical-inferences).", tags$br(), tags$br(),
                    "To prevent Lindley's paradox, one would need to lower the alpha level as a function of the statistical power. Good (1992) notes: 'we have empirical evidence that sensible P values are related to weights of evidence and, therefore, that P values are not entirely without merit. The real objection to P values is not that they usually are utter nonsense, but rather that they can be highly misleading, especially if the value of N is not also taken into account and is large.'
                    Therefore, we justify the alpha level as a function of sample size by relating it to Bayes Factors. A Bayes factor compares the likelihood of the data under the alternative hypothesis and under the null hypothesis. Therefore, setting the alpha level to always correspond to at least Bayes factor 1 avoids the Lindley paradox.
                    However, in Bayesian statistics, a Bayes factor of 1 or large is only regarded as weak evidence and we might wish to, for example, achieve at least moderate evidence if the p-value is significant. Therefore, we can adjust the desired evidence by using the slider."
                    )
              )
              ),
      tabItem(tabName = "anova",
            fluidRow(
              box(
                title = "Input parameters",
                numericInput("df1", "df1", 1, min = 0),
                numericInput("df2", "df2", 10, min = 0),
                sliderInput("evidence2", "How much more likely should the data at least be under the alternative hypothesis?",
                            min = 1,
                            max = 10,
                            value = 3),
                textOutput("likelyanova"),
                br(),
                selectInput("paired", "Within or Between Subjects?",
                            c("within" = TRUE,
                              "between" = FALSE
                            )),
                actionButton("power_start4", "Calculate")
              ),
            infoBoxOutput("anovabox"),
            box(plotOutput("plotanova")),
            box(title = "Explanation",
                status = "warning",
                solidHeader = TRUE, collapsible = TRUE,
                "The idea behind this recommendation is discussed most extensively by Leamer, 1978. He writes 'The rule of thumb quite popular now, that is, setting the significance level arbitrarily to .05, is shown to be deficient in the sense that from every reasonable viewpoint the significance level should be a decreasing function of sample size.' This was already recognized by Jeffreys (1939), who discusses ways to set the alpha level in Neyman-Pearson statistics: 'We should therefore get the best result, with any distribution of alpha, by some form that makes the ratio of the critical value to the standard error increase with n. It appears then that whatever the distribution may be, the use of a fixed P limit cannot be the one that will make the smallest number of mistakes.'", tags$br(), tags$br(),
                "The goal is to prevent Lindley's paradox (https://en.wikipedia.org/wiki/Lindley%27s_paradox). This is explained in more detail in week 1 of Daniel's MOOC (https://www.coursera.org/learn/statistical-inferences).", tags$br(), tags$br(),
                "To prevent Lindley's paradox, one would need to lower the alpha level as a function of the statistical power. Good (1992) notes: 'we have empirical evidence that sensible P values are related to weights of evidence and, therefore, that P values are not entirely without merit. The real objection to P values is not that they usually are utter nonsense, but rather that they can be highly misleading, especially if the value of N is not also taken into account and is large.'
                Therefore, we justify the alpha level as a function of sample size by relating it to Bayes Factors. A Bayes factor compares the likelihood of the data under the alternative hypothesis and under the null hypothesis. Therefore, setting the alpha level to always correspond to at least Bayes factor 1 avoids the Lindley paradox.
                However, in Bayesian statistics, a Bayes factor of 1 or large is only regarded as weak evidence and we might wish to, for example, achieve at least moderate evidence if the p-value is significant. Therefore, we can adjust the desired evidence by using the slider."
          )
        )
        ),
      # Third tab content
      tabItem(tabName = "about",
              h2("Justify Your Alpha: A Practical Guide"),
              h4("For an explanation why researchers should justify their alpha levels, see:"),
              h4("Lakens, D., Adolfi, F. G., Albers, C. J., Anvari, F., Apps, M. A. J., Argamon, S. E., … Zwaan, R. A. (2018). Justify your alpha. Nature Human Behaviour, 2, 168–171. https://doi.org/10.1038/s41562-018-0311-x"),
              h4("You can download the pre-print of this article at ", a("PsyArXiV", href="https://psyarxiv.com/9s3y6/")),
              h4("For a short introduction in why to lower your alpha level as a function of the sample size, see my ", a("blog post", href="http://daniellakens.blogspot.com/2018/12/testing-whether-observed-data-should.html"), ". For a short introduction on why and how to balance or minimize error rates, see my ", a("other blog post", href="http://daniellakens.blogspot.com/2019/05/justifying-your-alpha-by-minimizing-or.html"),"."),
              h4("Get the code at ", a("GitHub", href="https://github.com/Lakens/JustifieR")),
              h4("The best way to cite this app and the explanations of how to justify alpha levels in practice is through the preprint:"),
              h4("Maier & Lakens (2021). Justify Your Alpha: A Primer on Two Practical Approaches")
      )

    )
)
)



server <- function(input, output) {

  output$outH1H0 <- renderText({
    paste("This means that H1 is ", input$priorH1H0, " times as likley as H0 and that the probability of H1 is ", round(input$priorH1H0/(input$priorH1H0+1), digits = 2), sep = "")
  })

  output$outH1H02 <- renderText({
    paste("This means that H1 is ", input$priorH1H02, " times as likley as H0 and that the probability of H1 is ", round(input$priorH1H02/(input$priorH1H02+1), digits = 2), sep = "")
  })

  output$outT1T2 <- renderText({
    paste("This means that a false positive is ", input$costT1T2, " times as costly as a false negative.", sep = "")
  })

  output$outT1T22 <- renderText({
    paste("This means that a false positive is ", input$costT1T22, " times as costly as a false negative.", sep = "")
  })

  output$likelyttest <- renderText({
    if (input$evidence < 3)   {
      paste("This means that the data is at least ", input$evidence, " times more likely under the alternative than under the null. This avoids Lindleys paradox but implies only weak evidence.", sep = "")
    }
    else {
      if (input$evidence == 10){
        paste("This means that the data is at least ", input$evidence, " times more likely under the alternative than under the null. This corresponds to strong evidence for the alternative.", sep = "")
      } else{
        paste("This means that the data is at least ", input$evidence, " times more likely under the alternative than under the null. This corresponds to moderate evidence for the alternative.", sep = "")
      }
    }
  })

  output$likelyanova <- renderText({
    if (input$evidence2 < 3)   {
      paste("This means that the data is at least ", input$evidence2, " times more likely under the alternative than under the null. This avoids Lindleys paradox but implies only weak evidence.", sep = "")
    }
    else {
      if (input$evidence2 == 10){
        paste("This means that the data is at least ", input$evidence2, " times more likely under the alternative than under the null. This corresponds to strong evidence for the alternative.", sep = "")
      } else{
      paste("This means that the data is at least ", input$evidence2, " times more likely under the alternative than under the null. This corresponds to moderate evidence for the alternative.", sep = "")
      }
    }
  })


  observeEvent(input$power_start, {
    error <- isolate(input$error)
    power_function <- isolate(input$power_function)
    costT1T2 <- isolate(input$costT1T2)
    priorH1H0 <- isolate(input$priorH1H0)
    res <- optimal_alpha(power_function, costT1T2, priorH1H0, error, printplot = T)

    beta1 <- round(res$beta, digits = 4)
    alpha1 <- round(res$alpha, digits = 4)
    errorrate <- round(res$errorrate, digits = 4)

    costT1T22  <- input$costT1T22
    priorH1H02 <- input$priorH1H02

    # list(alpha1 = format(alpha1, digits = 10, nsmall = 5, scientific = FALSE),
    #      beta1 = format(beta1, digits = 10, nsmall = 5, scientific = FALSE))

    output$alpha1Box <- renderInfoBox({
      infoBox(
        "Alpha", alpha1,icon = icon("alpha"),
        color = "purple"
      )
    })

    output$beta1Box <- renderInfoBox({
      infoBox(
        "Beta", beta1, icon = icon("beta"),
        color = "green"
      )
    })

    output$error1Box <- renderInfoBox({
      infoBox(
        "Weighted Combined Error Rate", errorrate, icon = icon(""),
        color = "red"
      )
    })
    output$plot1 <- renderPlot({
      res$plot
    })

  })

  observeEvent(input$power_start2, {
    showModal(modalDialog("Estimating sample size, alpha level, and power. Please be patient, this might take several minutes.", footer=NULL))
    isolate(input$power_start2)
    error2 <- isolate(input$error2)
    errorrate2 <- isolate(input$errorrate2)
    power_function2 <- isolate(input$power_function2)
    costT1T22 <- isolate(input$costT1T22)
    priorH1H02 <- isolate(input$priorH1H02)
    res2 <- optimal_sample(power_function2, errorrate2, costT1T22, priorH1H02, error2, printplot = T)
    alpha2 <- round(res2$alpha, digits = 4)
    beta2 <- round(res2$beta, digits = 4)
    errorrates <- round(res2$errorrate, digits = 4)
    sample2 <- res2$samplesize

    output$alpha1Box2 <- renderInfoBox({
      infoBox(
        "Alpha", alpha2,icon = icon("alpha"),
        color = "purple"
      )
    })

    output$beta1Box2 <- renderInfoBox({
      infoBox(
        "Beta", beta2, icon = icon("beta"),
        color = "green"
      )
    })

    output$error1Box2 <- renderInfoBox({
      infoBox(
        "Weighted Combined Error Rate", errorrates, icon = icon(""),
        color = "red"
      )
    })

    output$sampleBox2 <- renderInfoBox({
      infoBox(
        "Sample Size", sample2, icon = icon(""),
        color = "yellow"
      )
    })

    output$plot2 <- renderPlot({
      res2$plot
    })
    removeModal()
  })
  # stats <- reactive({

  # })
  observeEvent(input$power_start3, {

    evidence <- as.numeric(isolate(input$evidence))
    n1 <- isolate(input$n1)
    n2 <- isolate(input$n2)
    one.sided <- as.logical(isolate(input$one.sided))
    rscale <- isolate(input$rscale)

  output$ttestbox <- renderInfoBox({

    infoBox(
      "Alpha", paste0(round(ttestEvidence(as.numeric(evidence), n1, n2, as.logical(one.sided), rscale = rscale)[[1]], digits = 3)),
      icon = icon("alpha"),
      color = "purple"
    )
  })

  output$plotttest <- renderPlot({
    #ttestEvidence(as.numeric(input$evidence), input$n1, input$n2, as.logical(input$one.sided), printplot = T)




    lindley  <- ttestEvidence(1,   n1, n2 = n2, one.sided, rscale = rscale, printplot =F)[[1]]
    moderate <- ttestEvidence(3,   n1, n2 = n2, one.sided, rscale = rscale, printplot =F)[[1]]
    strong   <- ttestEvidence(10,  n1, n2 = n2, one.sided, rscale = rscale, printplot =F)[[1]]
    indicated <- ttestEvidence(evidence,  n1, n2 = n2, one.sided, rscale = rscale, printplot =F)[[1]]

    loops <- seq(from = 0, to = 7, by = 0.01)
    p <- numeric(length(loops))
    bf <- numeric(length(loops))
    #d <- numeric(length(loops))
    tval <- numeric(length(loops))
    i <- 0
    for(t in loops){
      i <- i+1
      if(one.sided){
        bf[i] <- exp(BayesFactor::ttest.tstat(t, n1, n2, rscale = rscale, nullInterval = c(0, Inf))$bf)
        if(n2 != 0){
          p[i] <- pt(t, ((n1+n2) - 2), lower=FALSE)
        } else {
          p[i] <- pt(t, (n1 - 1), lower=FALSE)
        }
      } else {
        bf[i] <- exp(BayesFactor::ttest.tstat(t, n1, n2, rscale = rscale)$bf)
        if(n2 != 0){
          p[i] <- 2*pt(t, ((n1+n2) - 2), lower=FALSE)
        } else {
          p[i] <- 2*pt(t, (n1 - 1), lower=FALSE)
        }
      }
      tval[i] <- t
      #d[i] <- t * sqrt((1/n1)+(1/n2))
    }
    plot(p, bf, type="l", lty=1, lwd=3, xlim = c(0, max(0.05, lindley)), ylim = c(0.1, 10), axes = F, xlab = "p-value", ylab = "Bayes factor", log = "y")
    axis(side=1, at = c(0, as.numeric(lindley), as.numeric(moderate), as.numeric(strong), 0.05, indicated), labels = c(0, round(lindley, digits = 3), round(moderate, digits = 3), round(strong, digits = 3), 0.05, round(indicated, digits = 3)),  lwd = 3, las = 3)
    axis(side=2, at = c(0.1, 0.33, 1, 3, 10), labels = c("1/10", "1/3", 1, 3, 10), lwd = 3)
    points(indicated, evidence, col = "red", lwd = 4)
    abline(h = c(0.1, 0.33, 1, 3, 10), col = "gray", lty = 2)
    abline(v = c(lindley, moderate, strong), lty = 3)
    abline(v = indicated, lty = 3, col = "red")
  })
  })

  observeEvent(input$power_start4, {

    evidence <- as.numeric(isolate(input$evidence2))
    df1      <- isolate(input$df1)
    df2      <- isolate(input$df2)
    paired   <- isolate(input$paired)

  output$anovabox <- renderInfoBox({

    infoBox(
      "Alpha", paste0(round(ftestEvidence(evidence, df1, df2, paired)[[1]], digits = 3)),
      icon = icon("alpha"),
      color = "purple"
    )
  })

  output$plotanova <- renderPlot({
    lindley  <- ftestEvidence(1, df1, df2, paired, printplot = F)[[1]]
    moderate <- ftestEvidence(3, df1, df2, paired, printplot = F)[[1]]
    strong   <- ftestEvidence(10, df1, df2, paired, printplot = F)[[1]]
    indicated<- ftestEvidence(evidence, df1, df2, paired, printplot = F)[[1]]

    loops <- seq(from = 0, to = 100, by = 0.01)
    p <- numeric(length(loops))
    bf <- numeric(length(loops))
    #d <- numeric(length(loops))
    fval <- numeric(length(loops))
    i <- 0
    for(f in loops){
      i <- i+1
      bf[i] <- bf_bic(f, df1, df2, paired)
      p[i] <- (1 - pf(f, df1, df2))
      fval[i] <- f
      #d[i] <- t * sqrt((1/n1)+(1/n2))
    }
    plot(p, bf, type="l", lty=1, lwd=3, xlim = c(0, max(0.05, lindley)), ylim = c(0.1, 10), axes = F, xlab = "p-value", ylab = "Bayes factor", log = "y")
    axis(side=1, at = c(0, as.numeric(lindley), as.numeric(moderate), as.numeric(strong), 0.05, indicated), labels = c(0, round(lindley, digits = 3), round(moderate, digits = 3), round(strong, digits = 3), 0.05, round(indicated, digits = 3)),  lwd = 3, las = 3)
    axis(side=2, at = c(0.1, 0.33, 1, 3, 10), labels = c("1/10", "1/3", 1, 3, 10), lwd = 3)
    points(indicated, evidence, col = "red", lwd = 4)

    abline(h = c(0.1, 0.33, 1, 3, 10), col = "gray", lty = 2)
    abline(v = c(lindley, moderate, strong), lty = 3)
    abline(v = indicated, lty = 3, col = "red")
  })
  })

}


# Run the application
shinyApp(ui, server)
