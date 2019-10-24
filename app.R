## app.R ##
# library(tidyverse)
library(shiny)
library(shinydashboard)
library(stringr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(DT)
library(broom)
library(ggthemes)
library(readxl)
# library(ggnewscale)

# read in input datasets
lb_df_in <- read.delim("data/SIVcpzLB715_variant_table.txt", sep="\t", header=TRUE)
ek_df_in <- read.delim("data/SIVcpzEK505_variant_table.txt", sep="\t", header=TRUE)
mb_df_in <- read.delim("data/SIVcpzMB897_variant_table.txt", sep="\t", header=TRUE)
sm_df_in <- read.delim("data/SIVsm_variant_table.txt", sep="\t", header=TRUE)
mac_df_in <- read.delim("data/SIVmac239_variant_table.txt", sep="\t", header=TRUE)
hu_df_in <- read.delim("data/SIVhu_variant_table.txt", sep="\t", header=TRUE)
b670_df_in <- read.delim("data/SIVB670_variant_table.txt", sep="\t", header=TRUE)

# read in data related to frequencies of amino acids in HIV vs. SIV
ppm_df_in <- read.delim("data/all_ppm_data.txt", sep="\t", header=TRUE)
map_df_in <- read.delim("data/protein_msa_mapping.txt", sep="\t", header=TRUE)

# format ppm dataset
colnames(ppm_df_in) <- c("virus", "gene", "msa_position", "residue", "frequency") 
# replace gen "NA" values for stock with 0
ppm_df_in$frequency[is.na(ppm_df_in$frequency)] <- 0
ppm_df_in$frequency <- round(ppm_df_in$frequency, 2)
# format mapping dataset
colnames(map_df_in) <- c("siv_strain", "gene", "native_position", "msa_position")

ppm_df_spread <- spread(ppm_df_in, virus, frequency)

# join tables to create df by strain/position/amino acid & frequencies in HIV & SIV
ppm_df <- right_join(map_df_in, ppm_df_spread) 

# read in metadata.  For each dataset, this is the virus, the generation, the week, and the replicate (A or B)
metadata <- read_excel("data/metadata.xlsx")

# convert & tidy datasets
pre_process_dataset <- function(df){
  # get rid of water control datasets
  df <- df %>% select(-contains("Water"))
  
  df <- df  %>% gather(key = dataset_id, value = frequency, 7:ncol(df))
  
  # parse info out of dataset IDs
  # df$dataset_id <- str_replace_all(df$dataset_id, "Stock_virus", "stock_W0")
  
  # do a left join with metadata to 
  df <- left_join(df, metadata, by="dataset_id")
  
  # replace rep "envelope" values for cds with "env"
  df$cds <- str_replace(df$cds, "envelope", "env")
  
  # convert NA frequencies to 0s  
  df$frequency[is.na(df$frequency)] <- 0
  
  # pwg is position/week/generation ->  a unique identifier 
  df <- mutate(df, pwg = paste(as.character(position), week, gen))
  
  # variant name is position + aa change, e.g. 2832 gag K32N
  df <- mutate(df, variant_name = ifelse(is.na(cds), 
                                         position,
                                         paste0(position, " " , cds, " ", ref_aa, codon_number, alt_aa)))
  
  # variant name no position is e.g. gag K32N
  df <- mutate(df, variant_name_no_position = ifelse(is.na(cds), 
                                         position,
                                         paste0(cds, " ", ref_aa, codon_number, alt_aa)))
  
  
  # timepoints of sampling based on generation, week, and replicate (A or B)
  df <- mutate(df, timepoint = paste0("G", gen, "W", week, rep))
  
  # this allows plotting on a quasi-time x-axis, where generation 2 starts at week #30
  plot_week <- function(gen,week){
    weeks_per_gen <- 30
    # ifelse(gen==0, 0, 30+as.numeric(week), as.numeric(week))
    ifelse(gen==0, 0, (week + (weeks_per_gen * (as.integer(gen)-1))))
  }
  df$plot_week <- plot_week(df$gen, df$week)
  
  # now we need to create pseudo data points with plot_week = 30, and frequency = NA 
  # this will create a discontinuity in the plot to accomodate for the inter-generational
  # bottleneck...
  #
  # TODO: do this for all inter-generation transitions
  #
  add_pseudo_timepoint <- 0
  if (add_pseudo_timepoint) {
    # make pseudo timepoint rows
    df_stock <- filter(df, gen == 0)
    df_stock$plot_week[] <- 30
    df_stock$frequency[] <- NA
    df_stock$timepoint[] <- paste0("G1W30A")
    df <- rbind(df, df_stock)
    df_stock$rep[] <- "B"
    df_stock$timepoint[] <- paste0("G1W30B")
    df <- rbind(df, df_stock)
  }
  
  return (df)
}

# convert datasets to long (tidy) format  
lb_df_in <- pre_process_dataset(lb_df_in)
ek_df_in <- pre_process_dataset(ek_df_in)
mb_df_in <- pre_process_dataset(mb_df_in)
sm_df_in <- pre_process_dataset(sm_df_in)
mac_df_in <- pre_process_dataset(mac_df_in)
hu_df_in <- pre_process_dataset(hu_df_in)
b670_df_in <- pre_process_dataset(b670_df_in)
  
# merge frequency datasets with ppm data (add ppm data)

lb_df_in <- left_join(lb_df_in, filter(ppm_df, siv_strain == "LB715"), by = c("cds" = "gene", "codon_number" = "native_position", "ref_aa" = "residue"))  %>% 
  select (-siv_strain, -msa_position, ref_hiv=hiv, ref_siv=siv) 
lb_df_in <- left_join(lb_df_in, filter(ppm_df, siv_strain == "LB715"), by = c("cds" = "gene", "codon_number" = "native_position", "alt_aa" = "residue"))  %>% 
  select (-siv_strain, alt_hiv=hiv, alt_siv=siv) 
lb_df_in[]

ek_df_in <- left_join(ek_df_in, filter(ppm_df, siv_strain == "EK505"), by = c("cds" = "gene", "codon_number" = "native_position", "ref_aa" = "residue"))  %>% 
  select (-siv_strain, -msa_position, ref_hiv=hiv, ref_siv=siv) 
ek_df_in <- left_join(ek_df_in, filter(ppm_df, siv_strain == "EK505"), by = c("cds" = "gene", "codon_number" = "native_position", "alt_aa" = "residue"))  %>% 
  select (-siv_strain, alt_hiv=hiv, alt_siv=siv) 

mb_df_in <- left_join(mb_df_in, filter(ppm_df, siv_strain == "MB897"), by = c("cds" = "gene", "codon_number" = "native_position", "ref_aa" = "residue"))  %>% 
  select (-siv_strain, -msa_position, ref_hiv=hiv, ref_siv=siv) 
mb_df_in <- left_join(mb_df_in, filter(ppm_df, siv_strain == "MB897"), by = c("cds" = "gene", "codon_number" = "native_position", "alt_aa" = "residue"))  %>% 
  select (-siv_strain, alt_hiv=hiv, alt_siv=siv) 

# TODO: this for SIVsm, etc...


# setup UI and interactive inputs
ui <- dashboardPage(
    dashboardHeader(title = "Viral Variant Explorer"),
    dashboardSidebar(

    # option inputs for filtering
    
    h4("Filtering options:"),
    checkboxInput("in_CDS", "Must be in CDS", value=TRUE),
    checkboxInput("N_only", "Only non-synonymous", value=TRUE),
    checkboxInput("present_in_last_gen", "Present in last generation", value=FALSE),
    checkboxInput("positive_regression_slope", "Positive regression slope", value=FALSE),

    hr(),
    
    sliderInput("min_dataset_slider", label="Present in at least N datasets: ", min=1, max=10, value=4),
    sliderInput("min_max_spread", label="Min. difference between minimum and maximum frequencies: ", min=0, max=1, value=0),
    sliderInput("min_frequency_last_timepoint", label="Min. mean frequency last timepoint: ", min=0, max=1, value=0.5),
    sliderInput("min_alt_hiv", label="Min. alt. hiv frequency: ", min=0, max=1, value=0.0),
    
    hr(),
    
    # output # of variants pre & post filtering
    strong("Total variants:"),
    verbatimTextOutput("pre_filtering_variants"),
      
    strong("Filtered variants:"),
    verbatimTextOutput("post_filtering_variants"),
    
    # checkboxInput("no_table", "Don't show table", value=FALSE),
    
    h4("Plotting options:"),
    # checkboxInput("one_plot", "All variants on one scatter plot", value=FALSE),
    checkboxInput("circles_not_lines", "Circles instead of lines", value=TRUE),
    checkboxInput("plot_labels", "Variant labels", value=F),
    
    hr(),
    h4("Download plots or data table:"),
    downloadButton('download_scatter_plot', label = "Download scatter plots"),
    hr(),
    downloadButton('download_genome_plot',  label = "Download genome plot"),
    hr(),
    downloadButton('download_data_table',  label = "Download data table")
    
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
                   "MB897" = "mb",
                   "B670" = "b670",
                   "Hu" = "hu",
                   "Mac239" = "mac",
                   "SM" = "sm"), inline = TRUE)
      
                   # "All" = "all" ), inline = TRUE)
    ),
    fluidRow(plotOutput("var_plot")),
    fluidRow(plotOutput("genome_plot")),
    fluidRow(column(dataTableOutput("var_table"), width=12, style = "font-size:80%"))
    
  ) # end dashboardBody
)

server <- function(input, output) { 
  
  rv <- reactiveValues()
  
  # this reaction to switching datasets
  observeEvent(input$dataset_input, {
    rv$data_in <- switch(input$dataset_input,
                         sm = sm_df_in,
                         lb = lb_df_in,
                         ek = ek_df_in,
                         mb = mb_df_in,
                         hu = hu_df_in,
                         mac = mac_df_in,
                         b670 = b670_df_in,
                         all = rbind(sm_df_in, lb_df_in, ek_df_in, mb_df_in),
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
    input$min_frequency_last_timepoint
    input$min_alt_hiv
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
    
    n_var <-nrow(filtered_df %>% group_by(variant_name) %>% summarise())
    print (paste(n_var, "variants post CDS"))
    
    # filter out synonymous variants (or not)
    if (input$N_only){
      filtered_df <- filtered_df %>% filter(N_or_S == "N")
    }
    
    n_var <-nrow(filtered_df %>% group_by(variant_name) %>% summarise())
    print (paste(n_var, "variants post N/S"))
    
    # filter variants whose value in the last datapoint(s) is above a certain value    
    last_gen = max(filtered_df$gen)
    in_last_gen <- filter(filtered_df, as.numeric(gen) == as.numeric(last_gen))
    print (paste0("last gen: ", last_gen))
    last_week = max(as.numeric(in_last_gen$week))
    print (paste0("last week: ", last_week))
    sufficient_freq_in_last_week <- filter(in_last_gen, as.numeric(week) == as.numeric(last_week)) %>% 
                 group_by(variant_name) %>%
                 summarize(
                   mean_freq = mean(frequency)
                 ) %>%
                 filter(mean_freq > input$min_frequency_last_timepoint)
    filtered_df <- filtered_df %>% filter(variant_name %in% sufficient_freq_in_last_week$variant_name)


    # filter variants whose maximum value is a certain amount higher than the minimum value
    min_max_pass <- filtered_df %>% filter (!is.na(frequency)) %>%
                    group_by(variant_name) %>% summarize (min_freq = min(frequency), max_freq = max(frequency), max_min = max_freq - min_freq) %>%
                    filter (max_min > input$min_max_spread)
    filtered_df <- filtered_df %>% filter(variant_name %in% min_max_pass$variant_name)
                    
    n_var <-nrow(filtered_df %>% group_by(variant_name) %>% summarise())
    print (paste(n_var, "variants post min/max"))
    
    # filter based on presence in min # of datasets
    # keep only those present in more than N obervations
    in_n <- filter(filtered_df, frequency > 0 & !is.na(frequency)) %>% 
      group_by(variant_name) %>% 
      summarise(
        n = n()
      ) %>%
      filter(n >= input$min_dataset_slider)
    
    filtered_df <- filtered_df %>% filter(variant_name %in% in_n$variant_name)
    
    n_var <-nrow(filtered_df %>% group_by(variant_name) %>% summarise())
    print (paste(n_var, "variants post min datasets"))
    
    # filter based on presence in last generation
    if (input$present_in_last_gen){
      max_gen = max(filtered_df$gen)
      
      in_max_gen <- filter(filtered_df, as.numeric(gen) == as.numeric(max_gen)) %>% 
       filter(frequency > 0 & !is.na(frequency)) %>%
       group_by(variant_name, gen) %>% 
       summarise(
         n = n()
       ) %>%
       filter(n > 0)
      
    
       filtered_df <- filtered_df %>% filter(variant_name %in% in_max_gen$variant_name)
    }
    
    n_var <-nrow(filtered_df %>% group_by(variant_name) %>% summarise())
    print (paste(n_var, "variants post last gen"))
    
    # filter based on regression slope
    if (input$positive_regression_slope) {
      # calculate a frequency vs. week linear regression and extract the slope
      # filter out variants with a 0 or negative slope
      df_reg <- filtered_df %>% filter(!is.na(frequency)) %>% group_by(variant_name, rep) %>% select(variant_name, rep, plot_week, frequency)
      # lm y~x
      df_lm <- df_reg %>% do(fitWeek = lm(frequency ~ plot_week, data = .))
      dfHourCoef = tidy(df_lm, fitWeek)
    
      # filter those variants w/ slope > 0
      df_high_slope <- dfHourCoef %>% filter(term == "plot_week" & estimate > 0)
      filtered_df <- filtered_df %>% filter(variant_name %in% df_high_slope$variant_name)
    }
    n_var <-nrow(filtered_df %>% group_by(variant_name) %>% summarise())
    print (paste(n_var, "variants post pos. reg."))
    
    # filter variants based on difference between frequency of alternate allele in HIV & SIV 
    # compendium alignments
    if("alt_hiv" %in% colnames(filtered_df))
    {
      filtered_df <- filtered_df %>% 
        filter(alt_hiv >= input$min_alt_hiv)
    }
    
    n_var <-nrow(filtered_df %>% group_by(variant_name) %>% summarise())
    print (paste(n_var, "variants post min alt siv"))
    
    
    # keep track of # of post-filtering variants
    post_var <- nrow(filtered_df %>% group_by(variant_name) %>% summarise())
    output$post_filtering_variants <- renderText({post_var})
    print (paste(post_var, "rows post-filtering"))
    
    return (filtered_df)
  })
    
    # Render the data table
  
    output$var_table <- DT::renderDataTable(dt_function()) 
      
    
    prepare_dt <- function(){
      
      if("alt_hiv" %in% colnames(rv$data))
      {
        # This code collects the data into a table (tidy long -> wide format) and outputs in the order of timepoint collection
        data_table_unspread <- rv$data %>% filter(!is.na(frequency)) %>%
          # select (position, variant_name, gen, week, rep, frequency) %>% 
          select (position, variant_name, gen, week, rep, frequency, ref_hiv, ref_siv, alt_hiv, alt_siv) %>% 
          arrange(gen, as.numeric(week), rep)  %>% 
          mutate (ordered_timepoints = interaction(gen,week,rep), gen = NULL, week = NULL, rep=NULL) 
        
        ordered_timepoint_names <- unique(as.character(data_table_unspread$ordered_timepoints))
        data_table <- spread(data_table_unspread, ordered_timepoints, frequency) 
        data_table <- data_table[,c("variant_name", "ref_hiv", "ref_siv", "alt_hiv", "alt_siv", ordered_timepoint_names)]
        return (data_table)
      }
      else
      {
        # ref/alt frequencies in HIV/SIV sequences not available (yet) for this virus
        
        # This code collects the data into a table (tidy long -> wide format) and outputs in the order of timepoint collection
        data_table_unspread <- rv$data %>% filter(!is.na(frequency)) %>%
          # select (position, variant_name, gen, week, rep, frequency) %>% 
          select (position, variant_name, gen, week, rep, frequency) %>% 
          arrange(as.integer(gen), as.numeric(week), rep)  %>% 
          mutate (ordered_timepoints = interaction(gen,week,rep), gen = NULL, week = NULL, rep=NULL) 
        
        print ("here 2")
        # print (data_table_unspread[,c("ordered_timepoints", "position")])
        print (slice(data_table_unspread, 965:975))
        print ("here 3")
        
        ordered_timepoint_names <- unique(as.character(data_table_unspread$ordered_timepoints))
        data_table <- spread(data_table_unspread, ordered_timepoints, frequency) 
        data_table <- data_table[,c("variant_name", ordered_timepoint_names)]
        return (data_table)
      }
      
    }
      
    
    # output$var_table <- DT::renderDataTable({ 
    dt_function <- function() { 
      
      # this will make variant frequencies colored like a heatmap 
      breaks <- seq(0, 1, 0.1)
      blue_colors <- round(seq(255, 125, length.out = length(breaks) + 1), 0) %>%
      {paste0("rgb(" , . , "," , . , ",255)")}
      red_colors  <- round(seq(255, 125, length.out = length(breaks) + 1), 0) %>%
      {paste0("rgb(255," , . , "," , . , ")")}
      
      data_table <- prepare_dt()
      
      # create the data table object
      dat <- DT::datatable(data_table, 
                       options = list(
                         autoWidth = TRUE,
                         pageLength = 50,
                         scrollX = TRUE,
                         lengthMenu = c(10, 20, 50, 100, 200)
                         # columnDefs = list(list(width = '20px', targets = "_all" ))
                       )) %>%
        # the tail here omits first column (variant name) from coloring scheme
        formatStyle(tail(names(data_table), -5), backgroundColor = styleInterval(breaks, blue_colors)) %>%
        formatStyle(head(tail(names(data_table), -1),4), backgroundColor = styleInterval(breaks, red_colors))
        return(dat)
     }
    
    output$var_plot <- renderPlot(scatter_plot_function())
     scatter_plot_function = function(){
       # data_to_plot <- rv$data %>% arrange(as.numeric(as.character(position)))
       data_to_plot <- rv$data 
       
       # TODO: setup so you can highlight based on selections in table (or by gene, etc.)
       # highlight values selected (filtered) in DataTable
       # input$pct_id_table_rows_all is a vector index of the rv$data df 
       # that indicates which rows are actually selected in the DataTable
       # selected_rows <- input$pct_id_table_rows_all
       # selected_rows <- input$var_table_rows_selected
       # data_to_plot <- rv$data[selected_rows, ]
       
       # Facet wrap orders facet by the variable you specify, and it can be a bit tricky to change that order
       # see: 
       # https://kohske.wordpress.com/2010/12/29/faq-how-to-order-the-factor-variables-in-ggplot2/
       # for help
       data_to_plot$variant_name <- reorder(data_to_plot$variant_name, as.numeric(as.character(data_to_plot$position)))
       data_to_plot$variant_name_no_position <- reorder(data_to_plot$variant_name_no_position, as.numeric(as.character(data_to_plot$position)))
       
       gdtp <- data_to_plot %>% group_by(variant_name) %>% summarise (maxpos=max(position))
       print (gdtp)
       # ggplot(data=rv$data[!is.na(rv$data$frequency),], aes(x=plot_week, y= frequency, group=interaction(rep, position), color=rep)) +
       
       ggplot(data=data_to_plot, aes(x=plot_week, y= frequency, group=interaction(position, rep), color=rep, order=as.numeric(as.character(position)))) +
         geom_point(shape=20) +
         
         # TODO: replicates are confused...
         # TODO: discontinuos lines between generations
         
         geom_line(linetype=3) + 
         # geom_line(data=filter(data_to_plot, gen == 1)) + 
         # geom_line(data=filter(data_to_plot, gen == 2)) +
         
         # geom_smooth() +
         # theme(legend.position="none", strip.text = element_text(size=12) ) +  #We don't want a giant legend with each position #
         
         # theme_classic() +
         
         xlab("week post initial infection") +
         ylab("variant frequency") +
         # scale_y_log10() +
         # coord_fixed(ratio=20) + 
         facet_wrap(~variant_name_no_position) +
         # facet_wrap(~variant_name) +
         theme_tufte() + 
         theme(text = element_text(size=12, family="sans"),
               # legend.position="none",
               axis.line=element_line(),
               strip.text=element_text(size=10, family="sans", colour = "grey50")) +
         annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf)+
         annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf)+
         # scale_x_continuous( breaks=c(0,0.5,1)) + 
         scale_y_continuous(breaks=c(0,0.5,1), labels=c("0", "0.5", "1"))
       # scale_x_continuous(limits=c(10,35)) + scale_y_continuous(limits=c(0,400))
       # strip.background = element_blank()) 
       
     }
    
    output$genome_plot <- renderPlot(genome_plot_function())
    
    genome_plot_function = function(){
      
      # copy df & make a fake y axis variable that is a composite of the week plus the replicate
      data_to_plot <- rv$data %>% filter(!is.na(frequency)) %>% mutate(y_plot_week = plot_week + ifelse(rep == "A", 0, 2) )
      max_plot_y = max(data_to_plot$y_plot_week)
      variant_labels <- data_to_plot %>% group_by(variant_name) %>% summarize(name=variant_name_no_position[1], position=position[1], label_y=max_plot_y + 3)
      p <- ggplot(data=data_to_plot) 
      if (input$circles_not_lines) {
        p <- p + 
          geom_point(aes(y=y_plot_week, x=position, fill=frequency), shape=21, size=2.5, alpha=0.8, color="black", stroke=0.2)  +
          scale_fill_gradient(low="#FFDFDF", high="#FF0000", breaks=c(0,1)) 
      } else {
        p <- 
          p + geom_segment(data=data_to_plot, aes(x=position, xend=position, y=y_plot_week, yend=y_plot_week+1.5, color=frequency), size=1.5 ) +
          scale_color_gradient(low="#FFDFDF", high="#FF0000", breaks=c(0,1)) 
      } 
      p <- p + 
        scale_x_continuous(limits=c(0,10000)) +
        theme_classic() +
        theme(text = element_text(size=10, family="sans")) + 
        
        xlab ("position in genome (nt)") + 
        ylab ("week / replicate / passage") 
      
      if (input$plot_labels){
        p <- p + geom_text(data=variant_labels, aes(x=position, y=label_y, label=name))
      }
      
      p
      
    }
          
          # NS / S w/ different color scales - would need ggnewscale library to work but doesn't...
          # non-synonymous
          # geom_segment(data=filter(data_to_plot, N_or_S != "S"), aes(x=position, xend=position, y=y_plot_week, yend=y_plot_week+1.5, color=frequency), size=1.5 ) +
          # scale_color_gradient(low="#C7C7C7", high="#676767", breaks=c(0,1)) +
          # new_scale_color() +   # same as `new_scale("color")`
          # synonymous
          # geom_segment(data=filter(data_to_plot, N_or_S == "S"), aes(x=position, xend=position, y=y_plot_week, yend=y_plot_week+1.5, color=frequency), size=1.5 ) +
          # scale_color_gradient(low="#FFDFDF", high="#FF0000", breaks=c(0,1)) +
    
      #my @color_ranges =  ( [ 0xFFDFDF, 0xFF0000] , 
                            #[ 0xC7C7C7, 0x676767] );
    
    output$download_scatter_plot <- downloadHandler(
      filename = 'variant_scatter_plot.pdf',
      content = function(file) {
        device <- function(..., width, height) {
          grDevices::pdf(..., width = width, height = height)
        }
        ggsave(file, plot = scatter_plot_function(), device = "pdf", dpi=300, height=3, width=6, units="in")
      })
    
    output$download_genome_plot <- downloadHandler(
      filename = 'variant_genome_plot.pdf',
      content = function(file) {
        device <- function(..., width, height) {
          grDevices::pdf(..., width = width, height = height)
        }
        ggsave(file, plot = genome_plot_function(), device = "pdf", dpi=300, height=2.5, width=6, units="in")
      })
 
    output$download_data_table <- downloadHandler(
      filename = 'data_table.tsv',
      content = function(file) {
        write.table(prepare_dt(), file, sep='\t', col.names = NA)
      }
    )
 
} # end server function

shinyApp(ui, server)


