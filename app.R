## app.R ##
library(tidyverse)
library(shiny)
library(shinydashboard)
# library(stringr)
# library(dplyr)
# library(tidyr)
# library(ggplot2)
library(DT)
library(broom)

# read in input datasets
lb_df_in <- read.delim("data/SIVcpzLB715_variant_table.txt", sep="\t", header=TRUE)
ek_df_in <- read.delim("data/SIVcpzEK505_variant_table.txt", sep="\t", header=TRUE)
mb_df_in <- read.delim("data/SIVcpzMB897_variant_table.txt", sep="\t", header=TRUE)

# getrid of water control datasets
lb_df_in <- lb_df_in %>% select(-contains("Water"))
ek_df_in <- ek_df_in %>% select(-contains("Water"))
mb_df_in <- mb_df_in %>% select(-contains("Water"))

# convert & tidy datasets
pre_process_dataset <- function(df){
  df <- df  %>% gather(key = dataset_id, value = frequency, 7:ncol(df))
  
  # parse info out of dataset IDs
  df$dataset_id <- str_replace_all(df$dataset_id, "Stock_virus", "stock_W0")
  df$week   <- str_match(df$dataset_id, "_W([0-9]{1,2})_")[,2]
  df$gen    <- str_match(df$dataset_id, "_G([12])_")[,2]
  # replace gen "NA" values for stock with 0
  df$gen[is.na(df$gen)] <- 0
  
  mice           <- str_match(df$dataset_id, "G[12]_([J0-9]{4,5})_W[0-9]{1,2}_")[,2]
  
  # rename replicate A/B from mouse IDs 
  mice <- str_replace_all(mice, "J2626", "A")
  mice <- str_replace_all(mice, "J2627", "B")
  
  mice <- str_replace_all(mice, "J2778", "A")
  mice <- str_replace_all(mice, "J2779", "B")
  
  mice <- str_replace_all(mice, "J2663", "A")
  mice <- str_replace_all(mice, "J2666", "B")
  
  mice <- str_replace_all(mice, "2275", "A")
  mice <- str_replace_all(mice, "J2796", "B")
  
  mice <- str_replace_all(mice, "2139", "A")
  mice <- str_replace_all(mice, "2143", "B")
  
  mice <- str_replace_all(mice, "2263", "A")
  mice <- str_replace_all(mice, "2264", "B")
  
  df$rep <- mice
  
  # replace rep "NA" values for stock with A
  df$rep[is.na(df$rep)] <- "A"
  
  # remove observations w/ NA frequencies 
  # df <- filter(df, !is.na(frequency))
  
  # convert NA frequencies to 0s
  df$frequency[is.na(df$frequency)] <- 0
  
  # pwg is position/week/generation ->  a unique identifier 
  df <- mutate(df, pwg = paste(as.character(position), week, gen))

  # variant name is e.g. gag K32N
  df <- mutate(df, variant_name = ifelse(is.na(cds), 
                                         position,
                                         paste0(position, " " , cds, " ", ref_aa, codon_number, alt_aa)))
  
  
  
  # timepoints of sampling based on generation, week, and replicate (A or B)
  df <- mutate(df, timepoint = paste0("G", gen, "W", week, rep))
  

  # this allows plotting on a quasi-time x-axis, where generation 2 starts at week #30
  plot_week <- function(gen,week){
    ifelse(gen==2, 30+as.numeric(week), as.numeric(week))
  }
  df$plot_week <- plot_week(df$gen, df$week)
  
  return (df)
}
  
  
# convert datasets to long (tidy) format  
lb_df_in <- pre_process_dataset(lb_df_in)
ek_df_in <- pre_process_dataset(ek_df_in)
mb_df_in <- pre_process_dataset(mb_df_in)
  

# setup UI and interactive inputs
ui <- dashboardPage(
  dashboardHeader(title = "Viral Variant Explorer"),
  dashboardSidebar(
    # option inputs for filtering
    checkboxInput("in_CDS", "Must be in CDS", value=TRUE),
    checkboxInput("N_only", "Only non-synonymous", value=TRUE),
    checkboxInput("present_in_last_gen", "Present in last generation", value=FALSE),
    checkboxInput("positive_regression_slope", "Positive regression slope", value=FALSE),

    sliderInput("min_dataset_slider", label="Present in at least N datasets: ", min=1, max=10, value=4),
    sliderInput("min_max_spread", label="Min. difference between minimum and maximum frequencies: ", min=0, max=1, value=0),
    
    # output # of variants pre & post filtering
    strong("Total variants:"),
    verbatimTextOutput("pre_filtering_variants"),
      
    strong("Filtered variants:"),
    verbatimTextOutput("post_filtering_variants"),
    
    # checkboxInput("no_table", "Don't show table", value=FALSE),
    
    downloadButton('download_file', label = "Download plot")
    
    
  ), #Sidebar
  dashboardBody(
      # todo: generate dataset radio buttons programatically, e.g.
      # lapply(1:5, function(i) {
        # selectInput(paste0('a', i), paste0('SelectA', i),
                    # choices = sample(LETTERS, 5))
    fluidRow(
      radioButtons("dataset_input", "Dataset:",
                 c("LB715" = "lb",
                   "EK505" = "ek",
                   "MB897" = "mb"), inline = TRUE)
    ),
    fluidRow(plotOutput("var_plot")),
    fluidRow(column(dataTableOutput("var_table"), width=12, style = "font-size:80%"))
    
  ) # end dashboardBody
)

server <- function(input, output) { 
  
  rv <- reactiveValues()
  
  # this reaction to switching datasets
  observeEvent(input$dataset_input, {
    rv$data_in <- switch(input$dataset_input,
                         lb = lb_df_in,
                         ek = ek_df_in,
                         mb = mb_df_in,
                         lb_df_in)
  })
  

  # these inputs changing will triggger filtering
  # this observeEvent controls the reactivity --> when any input changes, filter variants appropriately
  observeEvent({
    input$in_CDS
    input$N_only 
    input$present_in_last_gen
    input$positive_regression_slope
    input$min_dataset_slider
    input$min_max_spread
    rv$data_in
    },{
      rv$data <- filter_variants(rv$data_in)
  })
  
  # a function to actually do the filtering...
  filter_variants <- isolate(function(df) {
    print ("filtering variants!")
    
    # keep track of # of total variants
    pre_var <-nrow(df %>% group_by(variant_name) %>% summarise())
    output$pre_filtering_variants <- renderText({pre_var})
    print (paste(pre_var, "rows pre-filtering"))
    
    # start w/ the un-filtered df
    filtered_df <- df %>% arrange(position)
    
    # filter out variants not in CDS (or not)
    if (input$in_CDS){
      filtered_df <- filtered_df %>% filter(!is.na(cds))
    }
    
    # filter out synonymous variants (or not)
    if (input$N_only){
      filtered_df <- filtered_df %>% filter(N_or_S == "N")
    }
    
    # filter variants whose maximum value is a certain amount higher than the minimum value
    min_max_pass <- filtered_df %>% 
                    group_by(variant_name) %>% summarize (min_freq = min(frequency), max_freq = max(frequency), max_min = max_freq - min_freq) %>%
                    filter (max_min > input$min_max_spread)
    filtered_df <- filtered_df %>% filter(variant_name %in% min_max_pass$variant_name)
                    
    
    # filter based on presence in min # of datasets
    # keep only those present in more than N obervations
    in_n <- filter(filtered_df, frequency > 0) %>% 
      group_by(variant_name) %>% 
      summarise(
        n = n()
      ) %>%
      filter(n > input$min_dataset_slider)
    
    filtered_df <- filtered_df %>% filter(variant_name %in% in_n$variant_name)
    
    # filter based on presence in last generation
    if (input$present_in_last_gen){
      max_gen = max(filtered_df$gen)
      
      in_max_gen <- filter(filtered_df, as.numeric(gen) == as.numeric(max_gen)) %>% 
       filter(frequency > 0) %>%
       group_by(variant_name, gen) %>% 
       summarise(
         n = n()
       ) %>%
       filter(n > 0)
    
       filtered_df <- filtered_df %>% filter(variant_name %in% in_max_gen$variant_name)
    }
    
    # filter based on regression slope
    if (input$positive_regression_slope) {
      # calculate a frequency vs. week linear regression and extract the slope
      # filter out variants with a 0 or negative slope
      df_reg <- filtered_df %>% group_by(variant_name, rep) %>% select(variant_name, rep, plot_week, frequency)
      # lm y~x
      df_lm <- df_reg %>% do(fitWeek = lm(frequency ~ plot_week, data = .))
      dfHourCoef = tidy(df_lm, fitWeek)
    
      # filter those variants w/ slope > 0
      df_high_slope <- dfHourCoef %>% filter(term == "plot_week" & estimate > 0)
      filtered_df <- filtered_df %>% filter(variant_name %in% df_high_slope$variant_name)
    }
    
    # keep track of # of post-filtering variants
    post_var <- nrow(filtered_df %>% group_by(variant_name) %>% summarise())
    output$post_filtering_variants <- renderText({post_var})
    print (paste(post_var, "rows post-filtering"))
    
    return (filtered_df)
  })
    
    # Render the data table
  
    output$var_table <- DT::renderDataTable({ 
      
      # this will make variant frequencies colored like a heatmap
      brks <- seq(0, 1, 0.1)
      clrs <- round(seq(255, 125, length.out = length(brks) + 1), 0) %>%
      {paste0("rgb(" , . , "," , . , ",255)")}
      
      # This code collects the data into a table (tidy long -> wide format) and outputs in the order of timepoint collection
      data_table_unspread <- rv$data %>% 
        select (position, variant_name, gen, week, rep, frequency) %>% 
        arrange(gen, as.numeric(week), rep)  %>% 
        mutate (ordered_timepoints = interaction(gen,week,rep), gen = NULL, week = NULL, rep=NULL) 
      
      ordered_timepoint_names <- unique(as.character(data_table_unspread$ordered_timepoints))
      data_table <- spread(data_table_unspread, ordered_timepoints, frequency) 
      data_table <- data_table[,c("variant_name", ordered_timepoint_names)]
      
      # create the data table object
      dat <- DT::datatable(data_table, 
                       options = list(
                         autoWidth = TRUE,
                         pageLength = 50,
                         lengthMenu = c(10, 20, 50, 100, 200)
                         # columnDefs = list(list(width = '20px', targets = "_all" ))
                       )) %>%
        formatStyle(tail(names(data_table), -6), backgroundColor = styleInterval(brks, clrs))
        return(dat)
     })
    
    output$var_plot <- renderPlot(var_plot_function())
    
    var_plot_function = function(){
      ggplot(data=rv$data[!is.na(rv$data$frequency),], aes(x=plot_week, y= frequency, group=interaction(rep, position), color=rep)) +
        geom_point() +
        geom_line() + 
        # geom_smooth() +
        # theme(legend.position="none", strip.text = element_text(size=12) ) +  #We don't want a giant legend with each position #
        # scale_y_log10() +
        # coord_fixed(ratio=20) + 
        facet_wrap(~variant_name) 
    }
    
    output$download_file <- downloadHandler(
      filename = 'variant_plot.pdf',
      content = function(file) {
        device <- function(..., width, height) {
          grDevices::pdf(..., width = width, height = height)
        }
        ggsave(file, plot = var_plot_function(), device = device)
      })
 
  } # end server function

shinyApp(ui, server)

