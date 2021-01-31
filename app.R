library(ggplot2)
library(rgeos)
library(rgdal)
library(raster)
library(maptools)
library(sf)
library(dplyr)
library(tidyverse)
library(zoo)
library(ggpubr)
library(egg)
library(RColorBrewer)
library(gridExtra)
library(bayesplot)
library(rstanarm)
library(readxl)
library(shiny)
library(shinythemes)
library(shinycssloaders)
library(shinyhelper)
library(shinyjs)
library(shinyWidgets)
library(rstan)

# Sys.setlocale("LC_TIME", "English")
delta_epsilon<-21
ymax_kft <- 0.012
colors_Dark<-brewer.pal(7,"Dark2")
colors_Spectral<-brewer.pal(7,"Spectral")
font_size = 20
font_size_title = 20
lwd = 1
pt_size = 1
right_margin=1

ui<-fluidPage(
    theme = shinytheme("spacelab"),
    
    tabsetPanel(
        tabPanel(title = strong("One serology point"),
                 sidebarLayout(
                     sidebarPanel(
                         fileInput('datafile', 'Choose Template file',buttonLabel = "Upload template",accept=c('xlsx', '.xlsx')),
                         includeMarkdown("help_upload_template.md"),
                         
                         sliderInput("input_n", "Time lag between seroconvert and serorevert",min = 30, max = 730, value = 161),
                         numericInput("input_sero", "unadjusted seroprevalence", value = 0.297,min = 0,max = 1,step = 0.001),
                         numericInput("kse", "sensitivity", value = 1,min = 0,max = 1,step = 0.001),
                         numericInput("ksp", "specificity", value = 0.975,min = 0,max = 1,step = 0.001),
                         numericInput("N", "serological survey sample size", value = 9514),
                         dateInput('input_date',label = 'Date of serology survey',value = "2020-7-21"),
                         numericInput("chains", "MCMC chains", value = 1,min = 0,max = 5,step = 1),
                         numericInput("niter", "MCMC niter", value = 10000,min = 1000,max = 50000,step = 1000),
                         actionButton("goButton", "Go!",class = "btn-primary"),
                         downloadButton('downloadData', 'Download')),
                     
                     mainPanel(conditionalPanel(condition = "input.goButton", withSpinner(plotOutput('plot'), type = 6))))),
        # mainPanel(withSpinner(plotOutput('plot'))))),
        tabPanel(title = strong("Multiple serology points"),
                 sidebarLayout(
                     sidebarPanel(
                         fileInput('datafile2', 'Choose Template file',accept=c('xlsx', '.xlsx')),
                         includeMarkdown("help_upload_template2.md"),
                         numericInput("p", "population size", value = 9050506),
                         numericInput("chains2", "MCMC chains", value = 1,min = 0,max = 5,step = 1),
                         numericInput("niter2", "MCMC niter", value = 10000,min = 1000,max = 50000,step = 1000),
                         actionButton("goButton2", "Go!",class = "btn-primary"),
                         downloadButton('downloadData2', 'Download')),
                     
                     mainPanel(conditionalPanel(condition = "input.goButton2", withSpinner(plotOutput('plot2'), type = 6)),
                               conditionalPanel(condition = "input.goButton2", withSpinner(plotOutput('plot3'), type = 6)))))
    ))

server<-function(input, output,session) {
    
    Nsize<-reactive(input$N)
    unadj_sero<-reactive(input$input_sero)
    kd<-reactive(1/input$input_n)
    t0<-reactive(as.Date(input$input_date))
    n2<-reactive(as.integer(as.Date(input$input_date)-as.Date("2020-01-01")+1))
    niter<-reactive(input$niter)
    chains<-reactive(input$chains)
    
    death<-reactive({
        if (is.null(input$datafile)) return(NULL)
        data1<-read_excel(input$datafile$datapath, sheet = "Epidemics")
        names(data1) <- c("date", "death")
        data1$death
    })
    
    IFR<-reactive({
        if (is.null(input$datafile))return(NULL)
        data2<-read_excel(input$datafile$datapath, sheet = "Severity-Mortality")
        names(data2) <- c("Age","IFR")
        data2$IFR
    })
    
    P0<-reactive({
        if (is.null(input$datafile))return(NULL)
        data3<-read_excel(input$datafile$datapath, sheet = "Population")
        names(data3) <- c("Age","population")
        data3$population
    })
    
    observeEvent( input$goButton, {
        
        posi_nega <- c(rep(1,round(unadj_sero()*Nsize())),rep(0,Nsize()-round(unadj_sero()*Nsize())))
        sero_data <-  list(N =Nsize() , y = posi_nega)
        
        fit <-stan(file="SeroAdjust.stan",
                   data = sero_data,
                   iter = niter(),
                   chains = chains(),
                   control = list(adapt_delta = 0.99))
        posx <- rstan::extract(fit)$p
        kf<-sum(IFR()*P0())/sum(P0())
        death2<-death()[1:n2()]
        t2<-seq(1,n2(),by=1)
        if (is.null(input$datafile))return(NULL)
        
        n<-length(death())
        t<-seq(1,n,by=1)
        q<-seq(1:length(posx))
        epsilon<-seroprevalence<-matrix(0,length(posx),n)
        
        for (i in 1:length(posx)) {
            q[i]<-(1-kf)/kf*exp(-kd()*n2())*sum(exp(kd()*t2)*death2)/(posx[i]*sum(P0()))+sum(death2)/sum(P0())
            epsilon[i,]<-(1-kf)/kf*cumsum(death())/q[i]/(sum(P0())-cumsum(death())/q[i])
            seroprevalence[i,]<-(1-kf)/kf*exp(-kd()*t)*cumsum(exp(kd()*t)*death()/q[i])/(sum(P0())-death()/q[i])
        }
        
        epsilon<-epsilon[,(delta_epsilon+1):n]
        
        data1 = data.frame(output = c(rep("Exposure", n-delta_epsilon), rep("Seroprevalence", n)),
                           t=c(as.Date(as.Date("2019-12-31")+1:length(epsilon[1,])) ,as.Date(as.Date("2019-12-31")+1:length(seroprevalence[1,]))),
                           median = c(apply(epsilon, 2, function(x) quantile(x, probs = 0.5)),
                                      apply(seroprevalence, 2, function(x) quantile(x, probs = 0.5))),
                           lower1 = c(apply(epsilon, 2, function(x) quantile(x, probs = 0.025)),
                                      apply(seroprevalence, 2, function(x) quantile(x, probs = 0.025))),
                           upper1 = c(apply(epsilon, 2, function(x) quantile(x, probs = 0.975)),
                                      apply(seroprevalence, 2, function(x) quantile(x, probs = 0.975))),
                           lower2 = c(apply(epsilon, 2, function(x) quantile(x, probs = 0.25)),
                                      apply(seroprevalence, 2, function(x) quantile(x, probs = 0.25))),
                           upper2 = c(apply(epsilon, 2, function(x) quantile(x, probs = 0.75)),
                                      apply(seroprevalence, 2, function(x) quantile(x, probs = 0.75))))
        data2 = data.frame( t=t0(), value=mean(posx), upper= quantile(posx, probs = 0.75), lower = quantile(posx, probs = 0.25))
        
        output$plot <- renderPlot({
            ymax_exposure<-ceiling(max(data1$upper2*10))/10
            p1<-ggplot(data1, aes(x=t, y = median, group = output, colour = output)) +
                geom_line(size=lwd) +  ggtitle("Exposure & Seroprevalence")+
                # geom_ribbon(aes(ymin=lower1, ymax=upper1, fill = output), alpha=0.2, colour = NA)+
                geom_ribbon(aes(ymin=lower2, ymax=upper2, fill = output), alpha=0.5, colour = NA)+
                scale_y_continuous(breaks =  seq(from=0, to=ymax_exposure, by=0.1),limit = c(0, ymax_exposure))+
                scale_x_date(breaks = as.Date(c("2020-01-01","2020-03-01", "2020-05-01", "2020-07-01","2020-09-01", "2020-11-01","2021-01-01")), labels=c("Jan","Mar", "May", "Jul","Sep","Nov","Jan"), limit = as.Date(c("2020-01-01","2021-01-31")))+
                geom_pointrange(data=data2, aes(x=t,y=value,ymin=lower, ymax=upper), inherit.aes = FALSE, shape = 21, size=pt_size, colour = "black", fill = colors_Dark[2])
            styled1 <- p1 +
                scale_fill_brewer(palette = "Dark2")+
                scale_colour_brewer(palette = "Dark2")+
                theme_minimal() +
                ylab(" ") +
                xlab("  ")+
                theme(
                    text = element_text(size=font_size),
                    plot.title = element_text(face = "bold", size = font_size_title),
                    legend.background = element_rect(fill = "white", size = 1.25, colour = "white"),
                    legend.justification = c(0, 1),
                    legend.position = c(0.2,0.9),
                    legend.title = element_blank(),
                    axis.ticks = element_line(colour = "grey50", size = 0.2),
                    panel.grid.major = element_line(colour = "grey50", size = 0.2),
                    panel.grid.minor = element_blank(),
                    plot.margin = margin(t=0, r=right_margin, b=0, l=0, "cm")
                )
            styled1
        })
        output$downloadData <- downloadHandler(
            filename = function() {
                paste('OneSeroModelResults-', Sys.Date(), '.csv', sep='')
            },
            content = function(con) {
                colnames(data1)<-c("Label","Date","median","95%lower","95%upper","50%lower","50%upper")
                write.csv(data1, con,row.names = FALSE)
            }
        )
    })
    
    niter2<-reactive(input$niter2)
    chains2<-reactive(input$chains2)
    
    p<-reactive(input$p)
    death3<-reactive({
        if (is.null(input$datafile2)) return(NULL)
        data3<-read_excel(input$datafile2$datapath, sheet = "Epidemics")
        names(data3) <- c("date","adjsero","death","adjsero_lower","adjsero_upper")
        data3$death
    })
    
    deathNEW<-reactive({
        if (is.null(input$datafile2)) return(NULL)
        deathNEW<-read_excel(input$datafile2$datapath, sheet = "Epidemics")
        names(deathNEW) <- c("date","adjsero","death","adjsero_lower","adjsero_upper")
        deathNEW
    })
    
    adjsero<-reactive({
        if (is.null(input$datafile2)) return(NULL)
        data3<-read_excel(input$datafile2$datapath, sheet = "Epidemics")
        names(data3) <- c("date","adjsero","death","adjsero_lower","adjsero_upper")
        data3$adjsero
    })
    
    observeEvent( input$goButton2, {
        
        death_end<-which(!is.na(death3()))[length(which(!is.na(death3())))]
        daily_death<-death3()[1:death_end]
        sero<-round(p()*adjsero()[!is.na(adjsero())])
        n_days<-length(daily_death)
        cumul_death<-cumsum(daily_death)
        n_days2<-length(sero)
        t<-seq(1,n_days,by=1)
        t2<-t[!is.na(adjsero())]
        data4<- list( n_days2 = n_days2 ,n_days=n_days ,sero = sero,daily_death =daily_death,t=t,t2=t2)
        
        fit2 <-stan(file="SeroAdjustMulti.stan",
                    data = data4,
                    iter = niter2(),
                    chains = chains2(),
                    control = list(adapt_delta = 0.99))
        beta <- rstan::extract(fit2)$beta
        gamma<- rstan::extract(fit2)$gamma
        sim<-length(beta)
        x<-epsilon<-kft<-matrix(0,sim,n_days)
        posterior <- as.matrix(fit2)
        
        output$plot3 <- renderPlot({
            # color_scheme_set("gray")
            
            plot_title <- ggtitle("Posterior distributions",
                                  "with medians and 95% intervals")
            mcmc_areas(posterior,
                       pars = c("beta", "gamma"),
                       prob = 0.95) + plot_title+
                scale_y_discrete(labels=c(
                    "gamma" =  expression(InfectionFatalityRatio(IFR)),
                    "beta" = expression(SeroreversionRate)
                ))+
                theme_minimal() +
                theme(
                    text = element_text(size=font_size),
                    plot.title = element_text(face = "bold", size = font_size_title),
                    legend.background = element_rect(fill = "white", size = 1.25, colour = "white"),
                    legend.justification = c(0, 1),
                    legend.position = c(0.2,0.9),
                    legend.title = element_blank(),
                    axis.ticks = element_line(colour = "grey50", size = 1),
                    plot.margin = margin(t=0, r=right_margin, b=0, l=0, "cm"),
                    panel.grid.major = element_line(colour = "grey50", size = 0.2),
                    panel.grid.minor = element_blank()
                )
        })
        
        for (i in 1:sim) {
            x[i,]<-rnbinom(rep(1,n_days), size= 100, mu=cumsum(exp(beta[i]*t)*(1-gamma[i])/gamma[i]*daily_death)/(exp(beta[i]*t)))/(p()-cumul_death)
            epsilon[i,]<-cumsum((1-gamma[i])/gamma[i]*daily_death)/(p()-cumul_death)
            kft[i,]<-gamma[i]
        }
        epsilon<-epsilon[,(delta_epsilon+1):n_days]
        data1 = data.frame(output = c(rep("Exposure", n_days-delta_epsilon), rep("Seroprevalence", n_days)),
                           t=c(as.Date(as.Date("2020-01-01")+1:length(epsilon[1,])),as.Date(as.Date("2020-01-01")+1:length(x[1,]))),
                           median = c(apply(epsilon, 2, function(x) quantile(x, probs = 0.5)),
                                      apply(x, 2, function(x) quantile(x, probs = 0.5))),
                           lower1 = c(apply(epsilon, 2, function(x) quantile(x, probs = 0.025)),
                                      apply(x, 2, function(x) quantile(x, probs = 0.025))),
                           upper1 = c(apply(epsilon, 2, function(x) quantile(x, probs = 0.975)),
                                      apply(x, 2, function(x) quantile(x, probs = 0.975))),
                           lower2 = c(apply(epsilon, 2, function(x) quantile(x, probs = 0.25)),
                                      apply(x, 2, function(x) quantile(x, probs = 0.25))),
                           upper2 = c(apply(epsilon, 2, function(x) quantile(x, probs = 0.75)),
                                      apply(x, 2, function(x) quantile(x, probs = 0.75))))
        data2= data.frame( t=as.Date(deathNEW()$date)[1:n_days][t2], value= deathNEW()$adjsero[t2], lower= deathNEW()$adjsero_lower[t2], upper = deathNEW()$adjsero_upper[t2])
        
        output$plot2 <- renderPlot({
            
            ymax_exposure<-ceiling(max(data1$upper2*10))/10
            
            p1<-ggplot(data1, aes(x=t, y = median, group = output, colour = output)) +
                geom_line(size=lwd) +  ggtitle("Exposure & Seroprevalence")+
                geom_ribbon(aes(ymin=lower1, ymax=upper1, fill = output), alpha=0.2, colour = NA)+
                geom_ribbon(aes(ymin=lower2, ymax=upper2, fill = output), alpha=0.5, colour = NA)+
                scale_y_continuous(breaks =  seq(from=0, to=ymax_exposure, by=0.1),limit = c(0, ymax_exposure))+
                scale_x_date(breaks = as.Date(c("2020-01-01","2020-03-01", "2020-05-01", "2020-07-01","2020-09-01", "2020-11-01","2021-01-01")), labels=c("Jan","Mar", "May", "Jul","Sep","Nov","Jan"), limit = as.Date(c("2020-01-01","2021-01-31")))+
                geom_pointrange(data=data2, aes(x=t,y=value,ymin=lower, ymax=upper), inherit.aes = FALSE, shape = 21, size=pt_size, colour = "black", fill = colors_Dark[2])
            styled2<- p1 +
                scale_fill_brewer(palette = "Dark2")+
                scale_colour_brewer(palette = "Dark2")+
                theme_minimal() +
                ylab(" ") +
                xlab("  ")+
                theme(
                    text = element_text(size=font_size),
                    plot.title = element_text(face = "bold", size = font_size_title),
                    legend.background = element_rect(fill = "white", size = 1.25, colour = "white"),
                    legend.justification = c(0, 1),
                    legend.position = c(0.2,0.9),
                    legend.title = element_blank(),
                    axis.ticks = element_line(colour = "grey50", size = 0.2),
                    panel.grid.major = element_line(colour = "grey50", size = 0.2),
                    panel.grid.minor = element_blank(),
                    plot.margin = margin(t=0, r=right_margin, b=0, l=0, "cm")
                )
            styled2
        })
        output$downloadData2 <- downloadHandler(
            filename = function() {
                paste('MultiSeroModelResults-', Sys.Date(), '.csv', sep='')
            },
            content = function(con) {
                colnames(data1)<-c("Label","Date","median","95%lower","95%upper","50%lower","50%upper")
                write.csv(data1, con,row.names = FALSE)
            }
        )
    })
    
    
}
shinyApp(ui = ui, server = server)

