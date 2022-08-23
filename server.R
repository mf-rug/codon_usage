cust_rang <- function(x) {
  out <- c()
  for (i in x) {
    out <- c(out, 255 * sin(0.5* pi * (i**0.25)))
  }
  round(out,0)
}
gcode <- c("F", "F", "L", "L", "S", "S", "S", "S", "Y", "Y",
           "*", "*", "C", "C", "*", "W", "L", "L", "L", "L", "P", "P", "P",
           "P", "H", "H", "Q", "Q", "R", "R", "R", "R", "I", "I", "I", "M",
           "T", "T", "T", "T", "N", "N", "K", "K", "S", "S", "R", "R", "V",
           "V", "V", "V", "A", "A", "A", "A", "D", "D", "E", "E", "G", "G",
           "G", "G")
gcode <- Biostrings::AMINO_ACID_CODE[gcode]
gcode[is.na(gcode)] <- 'STP' 
names(gcode) <- c("TTT", "TTC", "TTA", "TTG", "TCT", "TCC",
                  "TCA", "TCG", "TAT", "TAC", "TAA", "TAG", "TGT", "TGC", "TGA",
                  "TGG", "CTT", "CTC", "CTA", "CTG", "CCT", "CCC", "CCA", "CCG",
                  "CAT", "CAC", "CAA", "CAG", "CGT", "CGC", "CGA", "CGG", "ATT",
                  "ATC", "ATA", "ATG", "ACT", "ACC", "ACA", "ACG", "AAT", "AAC",
                  "AAA", "AAG", "AGT", "AGC", "AGA", "AGG", "GTT", "GTC", "GTA",
                  "GTG", "GCT", "GCC", "GCA", "GCG", "GAT", "GAC", "GAA", "GAG",
                  "GGT", "GGC", "GGA", "GGG")
nuc <- c("T", "C", "A", "G")
snuc = c(as.vector(sapply(nuc, rep, 16)),
         rep(as.vector(sapply(nuc, rep, 4)), 4),
         rep(nuc, 16))
codons <- apply(matrix(snuc, ncol = 3), 1, function(i) paste(i, collapse = ""))
nucdf <- data.frame(x = as.factor(rep(1L:64, 3)),
                    y = c(rep(0.2, 64), rep(1.8, 64), rep(2.7, 64)),
                    h = c(rep(2, 64), rep(1.2, 64), rep(0.6, 64)))
nucdf$nuc <- snuc
nucdf$nuc2 <- paste0(nucdf$nuc, '1')
nucdf$nuc3 <- paste0(nucdf$nuc, '2')
nucdf$aa <- rep(unname(gcode[codons]),3)

server <- function(input, output, session) {
  if (!all(file.exists('Codon_usage_of_all_refseq_species.csv'), file.exists('Codon_usage_of_important_refseq_species.csv'))) {
    cat(file=stderr(),'Processed data base not found, creating files...')
    # read the file with codon usage from all species, downloaded from https://dnahive.fda.gov/dna.cgi?cmd=codon_usage&id=537&mode=cocoputs
    df <- read_tsv('o537-Refseq_species.tsv', show_col_types = FALSE) %>%as.data.frame()
    
    # select just the Taxid, Species, and Codon columns
    df <- df[,c(3:4, 13:ncol(df))]
    
    # calculate the mean across various assemblies with slightly changing values
    df <- df %>% group_by(Species) %>% mutate(across(TTT:GGG, mean))
    
    # keep just a single row for each individual species
    df <- df[!duplicated(df),]
    
    # write output
    write_csv(df, 'Codon_usage_of_all_refseq_species.csv', col_names = TRUE)
    cat(file=stderr(),'Wrote all species')
    write_csv(df[df$Taxid %in% sort(83333, 9606, 4922, 4932, 10090, 10117,3702, 6239, 7227, 7955,8355),],
              'Codon_usage_of_important_refseq_species.csv', col_names = TRUE)
    cat(file=stderr(),'Wrote important species')
    df_imp <- df
  } else {
    df_imp <- read_csv('Codon_usage_of_important_refseq_species.csv', col_names = TRUE, show_col_types = FALSE)
  }
  df_imp <- df_imp[order(df_imp$Species),]
  sp_id <- paste0(df_imp$Species, ' [Taxid:', df_imp$Taxid, ']')
  codon_to_aa <- Biostrings::GENETIC_CODE %>% as.data.frame()
  codon_to_aa$... <- Biostrings::AMINO_ACID_CODE[codon_to_aa$.]
  codon_to_aa[is.na(codon_to_aa)] <- '*'
  updateSelectizeInput(session, 'species', choices = levels(factor(sp_id)), selected = "Escherichia coli K-12 [Taxid:83333]", server = TRUE)

  df <- eventReactive(eventExpr = input$all_species, valueExpr = {
    if (input$all_species) {
      if (!exists("df_full")) {
        df_full <- read_csv('Codon_usage_of_all_refseq_species.csv', col_names = TRUE, show_col_types = FALSE)
      }
      df <- df_full
    } else {
      df <- df_imp
    }
    df[order(df$Species),]
  })
  
  observeEvent(input$all_species, {
    df <- df()
    sp_id <- paste0(df$Species, ' [Taxid:', df$Taxid, ']')
    updateSelectizeInput(session, 'species', choices = sp_id, selected = NA, server = TRUE)
  }, ignoreInit = TRUE)
  
  
  df_out <- eventReactive(input$species, {
    df <- df()
    taxid <- input$species %>% str_extract('[0-9][0-9]*(?=\\]$)')
    df_out <- df[df$Taxid == taxid,]
    df_out <- pivot_longer(df_out, cols = TTT:GGG)
    df_out$aa <- codon_to_aa[df_out$name, '...']
    df_out$a <- codon_to_aa[df_out$name, '.']
    df_out <- df_out %>% group_by(aa) %>% mutate(value_norm = value / sum(value))
    df_out
  })
  
  output$table <- renderDT({
    if (input$species != '') {
      df_out <- df_out()[,c(3,5,6,4,7)]
      df_out$value <- round(df_out$value, 0)
      df_out$value_norm <- round(df_out$value_norm, 3)
      colnames(df_out) <- c('Codon', 'Amino acid 3', 'Amino acid', 'Total count', 'Fraction')
      
      brks <- sort(df_out$Fraction)
      vals <- cust_rang(seq(min(df_out$Fraction), max(df_out$Fraction), length.out = length(brks) + 1))
      clrs <- paste0("rgb(255,", vals, ",", vals, ")")
      
      datatable(df_out,
                rownames= FALSE,
                options = list(pageLength = 100,
                               autoWidth = TRUE,
                               columnDefs = list(
                                 list(width = '6vw', targets = c("Codon", "Amino acid", "Amino acid 3")),
                                 list(width = '24vw', targets = c("Fraction"))
                               ),
                               dom = 'ft')) %>%
        formatStyle(0, target= 'row', lineHeight='70%') %>% 
        formatStyle(columns = 'Total count', 
                    background = styleColorBar(range(df_out$`Total count`), 'lightblue'),
                    backgroundSize = '98% 88%',
                    backgroundRepeat = 'no-repeat',
                    backgroundPosition = 'center') %>% 
        formatStyle('Fraction', backgroundColor = styleInterval(brks, clrs))
    }
  })
  
  output$plot <- renderPlot({
    df_out <- df_out()
    if (input$species != '') {
      p <- ggplot(df_out)
      if (input$norm) {
        p <- p + 
          aes(x= reorder(name, aa, sort), y = value_norm) +
          scale_y_continuous(limits = c(0, max(df_out$value_norm)), expand = c(0,0)) +
          ylab('Fraction')
      } else {
        p <- p + 
          aes(x= name, y = value) +
          scale_y_continuous(limits = c(0, max(df_out$value)), expand = c(0,0)) +
          ylab('Total count')
      }
      p <- p +
        geom_col(aes(fill =aa), 
                 show.legend = FALSE) +
        xlab('Codon') +
        # scale_fill_manual(values = unname(pals::polychrome())) +
        theme_bw() +
        theme(axis.title = element_text(size = s_f_d() * 19, face = "bold"),
              axis.text = element_text(size = s_f_d() * 16),
              axis.text.x = element_text(size = s_f_d() * 13, angle = 45, vjust = .85, hjust = .85),
              strip.text = element_text(size = s_f_d() * 14.5))
      
      if (input$by_aa) {
        p <- p + 
          facet_grid(~aa, scales = 'free_x', space = "free_x")
      } else {
        p <- p + geom_text(aes(label = aa),position = position_stack(vjust = 0.5))
      }
      p
    }
  })
  
  # get screen resolution to adjust font sizes by a factor
  s_f <- eventReactive(input$dimension,{
    input$dimension[1] / 2500
  })
  
  # prevent graph from being redrawn many times if user resizes window
  s_f_d <- debounce(s_f, 1000)
  
  output$plotsun <- renderPlot({
    if (!identical(s_f_d(), numeric(0))) {
  
      ggplot(nucdf, aes(x = x, y = y)) +
        # 3' big black points in the background
        geom_point(data = NULL, x = 0.5, y= 3.4, size = s_f_d() * 20, color = 'black', shape = 16) + 
        geom_point(data = NULL, x = 48.5, y= 3.4, size = s_f_d() * 20, color = 'black', shape = 16) +
        geom_point(data = NULL, x = 96.5, y= 3.4, size = s_f_d() * 20, color = 'black', shape = 16) +
        geom_point(data = NULL, x = 144.5, y= 3.4, size = s_f_d() * 20, color = 'black', shape = 16) +
        
        # AA tiles are the outermost ring
        geom_tile(data = nucdf[129:192,], aes(y = y + 0.24, height = h + 0.3, fill = aa), color = 'transparent', size = 0, show.legend = FALSE) +
        geom_vline(xintercept = c(2,4,8,10,12,14,15,16,20,24,26,28,32,35,36,40,42,44,46,48,52,56,58,60,64,66,68,72,74,76,78,79,80,84,88,90,92,96,99,100,104,106,108,110,
                                  112,116,120,122,124,128,130,132,136,138,140,142,143,144,148,152,154,156,160,163,164,168,170,172,174,176,180,184,186,188,192) + 0.5,
                   linetype = 1, size = 0.3, color = "black") +
        
        # first Nuc tiles are the outer ring
        geom_tile(data = nucdf[129:192,], aes(height = h, fill = nuc3), color = 'black', size = 0.5, show.legend = FALSE) +
        geom_vline(xintercept = seq(4,64, 4) + 0.5, linetype = 1, size = 1.1, color = "black") +
        
        
        # second tiles are the middle ring
        geom_tile(data = nucdf[65:128,], aes(height = h, fill = nuc2, color = nuc2), size = 1,  show.legend = FALSE) +
        geom_vline(xintercept = seq(4,66, 4) + 0.5,linetype = 1, size = 1.1, color = "black") +
        
        # third tiles are the inner ring
        geom_tile(data = nucdf[1:64,], aes(height = h, fill = nuc, color = nuc), size = 1, show.legend = FALSE) +
        geom_vline(xintercept = c(0.5, 96.5, 48.5, 144.5), linetype = 1, size = 1.4, color = "black") +
        
        # thin line to separate same amino acids across the entire chart
        geom_vline(xintercept = c(2,4,8,10,12,14,15,16,20,24,26,28,32,35,36,40,42,44,46,48,52,56,58,60,64,66,68,72,74,76,78,79,80,84,88,90,92,96,99,100,104,106,108,110,
                                  112,116,120,122,124,128,130,132,136,138,140,142,143,144,148,152,154,156,160,163,164,168,170,172,174,176,180,184,186,188,192) + 0.5, 
                   linetype = 1, size = 0.3, color = "black", alpha =0.06) +
        
        # border around amino acid ring
        geom_hline(yintercept = c(3,3.4), size = 1, color = 'black') +
        
        # amino acid radial text
        geom_textvline(aes(xintercept = as.numeric(x)), label = rep(unname(gcode[codons]),3), hjust = 0.975 + (s_f_d() * 0.02),
                       linetype = 0, color = "black", family = 'mono', fontface = 'bold', size = 4 * s_f_d()) +
        
        # 5' big black point, and label 5' and 3' points
        geom_point(data = NULL, x = 0, y= -1, size = s_f_d() * 20, color = 'black') +
        geom_text(data = NULL, aes(label = "3'"), x = 0.5, y= 3.4, size = s_f_d() * 6, color = 'white', vjust = -0.3) + 
        geom_text(data = NULL, aes(label = "3'"), x = 48.5, y= 3.4, size = s_f_d() * 6, color = 'white', hjust = 1.3) + 
        geom_text(data = NULL, aes(label = "3'"), x = 96.5, y= 3.4, size = s_f_d() * 6, color = 'white', vjust = 1.3) + 
        geom_text(data = NULL, aes(label = "3'"), x = 144.5, y= 3.4, size = s_f_d() * 6, color = 'white', hjust = -0.3) + 
        geom_text(data = NULL, x = 0, y= -1, label = "5'",  size = s_f_d() * 8, color = 'white') +
        
        # labels for the nuc rings
        geom_text(data = nucdf[129:192,], mapping = aes(label = nuc), size = s_f_d() * 6, fontface = 'bold') +
        geom_text(data = nucdf[seq(66, 126,length.out = 16),], mapping = aes(y= y, label = nuc), size = s_f_d() * 15, fontface = 'bold') +
        geom_text(data = nucdf[c(8,24,40,56),], mapping = aes(y= y, label = nuc), size = s_f_d() * 22, fontface = 'bold') +
        
        # axis
        coord_polar() +
        scale_y_continuous(limits = c(-1, 3.4)) +
        
        # colors
        scale_fill_manual(values = c("A" = input$g, "T" = input$h, "C" = input$i, "G" = input$j,
                                     "A1" = '#B97FFF', "T1" = '#FF75BC', "C1" = '#6BA9FF', "G1" = '#90FF8A',
                                     "A2" = '#D3ADFF', "T2" = '#FFADD6', "C2" = '#9CC5FF', "G2" = '#C0FFBD',
                                     'Phe'="grey55",'Leu'="grey70",'Ser'="green",'Tyr'="#ADBDAC",'STP'="white",
                                     'Cys'="#F4FC9F",'Trp'="grey50",'Pro'="grey70",'His'="#7759EB",'Gln'="#BEFFB2",
                                     'Arg'="blue",'Ile'="grey70",'Met'="#D4D6B4",'Thr'="#81FF6B",'Asn'="#95FF82",
                                     'Lys'="#5F57FF",'Val'="grey80",'Ala'="grey88",'Asp'="#FF5454",'Glu'="red",'Gly'="grey95")) +
        
        scale_color_manual(values = c("A" = input$g, "T" = input$h, "C" = input$i, "G" = input$j,
                                      "A1" = '#B97FFF', "T1" = '#FF75BC', "C1" = '#6BA9FF', "G1" = '#90FF8A',
                                      "A2" = '#D3ADFF', "T2" = '#FFADD6', "C2" = '#9CC5FF', "G2" = '#C0FFBD')) +
        theme(axis.ticks = element_blank(), axis.text = element_blank(),
              axis.title.y = element_blank(), axis.title.x = element_blank(),
              panel.border = element_blank(), panel.grid.major = element_blank(),
              rect = element_rect(fill = "transparent"),
              panel.background = element_rect(fill = "transparent",colour = NA),
              plot.background = element_rect(fill = "transparent",colour = NA),
              plot.margin = margin(s_f_d() * -1.8,s_f_d() * -1.8,s_f_d() * -1.8,s_f_d() * -1.8, "cm"))
    }
  }, bg="transparent")
  
}
