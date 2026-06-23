library(tidyverse)
library(karyoploteR)


##  This first section is to create individual sample ABH files from the main ABH file
##  Only run this section once to save time later. This section writes new files from the ABH file that are quicker to read later on.

# read in the ABH file (very large file that will take some time to read in)
abh_in = read.csv("input_data/tutorial_samples.ABH.csv")

# extract a vector of the sample names listed in the "id" column of the ABH file, excluding NA values
sample_names = abh_in$id[!is.na(abh_in$id)]

# create directory to store new individual sample ABH files in
dir.create(file.path("input_data/ABH_Individuals"))


# for loop to extract ABH file for each sample by its sample ID found in sample_names
for(sample_id in sample_names){
  abh_single_sample = filter(abh_in, grepl(sample_id, id)) %>% # filter to have only the current sample of interest
    select(id, starts_with("CM")) %>% # select the id column and all columns with ABH values by SNP
    pivot_longer(starts_with("CM"), names_to = "genome_location", values_to = "ABH") %>% # convert rows to columns
    separate(col = genome_location, into = c("chr_name", "chr_number", "position"), sep = "[.]") %>% # separate genome_location column to multiple columns
    unite(c("chr_name", "chr_number"), col="chr", sep=".") %>% # recombine chromosome name and number columns to one column THAT MUST BE NAMED "chr"
    mutate(position = as.numeric(position)) # save position of SNP as a numeric variable, rather than character
  
  # write single sample in ABH form as a single .csv file
  write.csv(abh_single_sample, paste0("input_data/ABH_Individuals/", sample_id, ".ABH.csv"), row.names = FALSE)
}


# save sample names for quick access later
write.csv(data.frame(sample_names), "input_data/ABH_Individual_sample_names.csv", row.names = F)

## << end of first section >>









## This section to the end of the code file is the process for making the introgression visualizations
# First, open the file with chromosome sizes for this species's genome. This should be downloadable online if the reference genome for your species of interest has been created
# Be sure that the path in read.table is correct for your reference genome's length file local storage location
genome_lengths_in = read.table("input_data/GCA_003086295.3_arahy.Tifrunner.gnm2.J5K5_genomic.length", header = F)

# add in a start position column, starting at 1 for each chromosome, for each chromosome. Also, rename the columns. THE COLUMNS MUST BE NAMED "chr", start", AND "end"
genome = genome_lengths_in %>%
  mutate(start=1, .before = V2) %>% 
  rename_with(~ c("chr", "end"), starts_with("V"))

# read in the sample names file to use in the upcoming for loop
sample_names = read.csv("input_data/ABH_Individual_sample_names.csv")# %>% 

# create a data frame from the genome lengths file with the start and end position of each chromosome. This will be used in later steps with the individual ABH files in the for loop
chr_start_end = genome %>% 
  pivot_longer(c(start,end), names_to="type", values_to = "position") %>% 
  mutate(ABH="T") %>% 
  select(-type)

# create a directory to store the final introgression plots in. (We will save them as jpeg images)
dir.create(file.path("Introgression Plots"))

# this for loop opens an individual sample ABH file one-by-one and creates an introgression plot that individual sample
for(sample in sample_names$sample_name){
  
  # indicate the beginning of the individual sample's plotting process
  print(paste(sample,"sample introgression plotting process started..."))
  
  # karyoploteR requires a data frame that indicates beginning and ending positions of sections in the chromosome that are A, B, or H. These next couple of blocks accomplish that
  sample_individual = read.csv(paste0("input_data/ABH_Individuals/", sample,".ABH.csv")) %>% 
    rbind(chr_start_end %>% mutate(id=sample, .before=chr)) %>% # combine individual sample's ABH file with the earlier data to also indicate where each chromosome begins and ends
    arrange(chr, position) %>% # now arrange all positions, including beginning and end of chromosomes, in order
    replace_na(list(ABH="N")) %>% # substitute NA values with an "N" character value
    filter(!str_detect(ABH, "N")) %>% # comment this out when wanting the NA snps included! Otherwise, this just ignores NA/"N" values
    mutate(ABH_prev = lag(ABH, n=1), ABH_next = lead(ABH, n=1)) %>% # create new columns of same ABH data, but offset to show next and previous SNP's ABH value
    mutate(is_start = case_when( # this switch case shows a position as a start position for an A, B, or H section depending on the previous SNP's value
      str_detect(ABH, "B") & str_detect(ABH_prev, "A|H|N|T") ~ "startB",
      str_detect(ABH, "H") & str_detect(ABH_prev, "A|B|N|T") ~ "startH",
      str_detect(ABH, "A") & str_detect(ABH_prev, "B|H|N|T") ~ "startA",
      str_detect(ABH, "N") & str_detect(ABH_prev, "B|H|A|T") ~ "startN",
      TRUE ~ "none" # in case previous SNP is 'T' (beginning of a chromosome)
    )) %>% 
    mutate(is_end = case_when( # this switch case shows a position as an end position for an A, B, or H sectioning depending on the value of the next SNP's value
      str_detect(ABH, "B") & str_detect(ABH_next, "A|H|N|T") ~ "endB",
      str_detect(ABH, "H") & str_detect(ABH_next, "A|B|N|T") ~ "endH",
      str_detect(ABH, "A") & str_detect(ABH_next, "B|H|N|T") ~ "endA",
      str_detect(ABH, "N") & str_detect(ABH_next, "B|H|A|T") ~ "endN",
      TRUE ~ "none" # in case next SNP is 'T' (end of a chromosome)
    )) %>% 
    select(-ABH_next, -ABH_prev) %>% # drop the offset columns now. Their purpose is done
    filter(str_detect(is_start, "start") | str_detect(is_end,"end")) # now drop all positions that are not the start or end of an A, B, or H section


    
  # determine start and stop SNPs of sections homozygous for A parent sections
  start_pos_a = sample_individual %>% 
    filter(is_start=="startA") %>% # filter to have only positions that are the start of an A section
    select(chr, position) %>%  # select only the chromosome and position of the start of an A section, dropping the is_start column (no longer needed and wanted)
    dplyr::rename(start = position) # rename the position to start. Name MUST BE "start"
  end_pos_a = sample_individual %>% # do same for end positions of A sections
    filter(is_end=="endA") %>% 
    select(position) %>% 
    dplyr::rename(end=position)
  
  # and for homozygous B parent sections
  start_pos_b = sample_individual %>% 
    filter(is_start=="startB") %>% 
    select(chr, position) %>% 
    dplyr::rename(start = position)
  end_pos_b = sample_individual %>% 
    filter(is_end=="endB") %>% 
    select(position) %>% 
    dplyr::rename(end=position)
  
  # and for sections heterozygous for both
  start_pos_h = sample_individual %>% 
    filter(is_start=="startH") %>% 
    select(chr, position) %>% 
    dplyr::rename(start = position)
  end_pos_h = sample_individual %>% 
    filter(is_end=="endH") %>% 
    select(position) %>% 
    dplyr::rename(end=position)
  
  # and for NA sections, if NAs have not been ignored
  start_pos_n = sample_individual %>% 
    filter(is_start=="startN") %>% 
    select(chr, position) %>% 
    dplyr::rename(start = position)
  end_pos_n = sample_individual %>% 
    filter(is_end=="endN") %>% 
    select(position) %>% 
    dplyr::rename(end=position)
  
  
  # combine dataframes with start and end position columns for each type of introgression (aka section of A, B, or H)
  a_parent_sections = cbind(start_pos_a, end_pos_a) %>% 
    unite(name, c(chr, start, end), remove=F) %>% 
    relocate(name, .after=end) %>% # name is needed for the format required by karyoploteR, but it doesn't do anything for how we are using karyoploteR
    mutate(gieStain="stalk") # the gieStain variable is needed for karyoploteR. The value of gieStain will correspond with a color
  b_parent_sections = cbind(start_pos_b, end_pos_b) %>% 
    unite(name, c(chr, start, end), remove=F) %>% 
    relocate(name, .after=end) %>% 
    mutate(gieStain="acen")
  h_sections = cbind(start_pos_h, end_pos_h) %>% 
    unite(name, c(chr, start, end), remove=F) %>% 
    relocate(name, .after=end) %>% 
    mutate(gieStain="acen2")
  na_sections = cbind(start_pos_n, end_pos_n) %>% # if NAs have been ignored, this dataframe will really be empty
    unite(name, c(chr, start, end), remove=F) %>% 
    relocate(name, .after=end) %>% 
    mutate(gieStain="gpos50")
  
  
  # based on our dataframe made from the genome LENGTH file, this will create sections for chromosome background. This dataframe file is in the same format as the A, B, H, and NA section dataframes above (because we are about to combine them all into one data frame)
  main_chrom_underlay = genome %>% 
    unite(name, c(chr, start, end), sep="-", remove=F) %>% 
    relocate(name, .after=end) %>% 
    mutate(gieStain="gneg")
  
  # here we combine all A, B, H, NA, and chromosomal underlay section data frames into one dataframe
  full_map = rbind(main_chrom_underlay, na_sections, a_parent_sections, b_parent_sections, h_sections) # to finish making a file like https://raw.githubusercontent.com/bernatgel/karyoploter_examples/master/Examples/Tutorial/CustomGenomes/mycytobands.txt

  
  

    
  
  # finally, the actual plotting
  
  # convert the genome file to one usable by karyoploteR
  peanut.genome = toGRanges(genome)
  
  # set some parameters for the plot
  pp = getDefaultPlotParams(plot.type=2)
  pp$ideogramheight = 400
  pp$leftmargin = .12

  # assemble the plot by layers
  kp = plotKaryotype(genome=peanut.genome, cytobands = full_map, chromosomes = "all", plot.params=pp, ideogram.plotter = NULL) # create the layout for the plot. It is sized based on peanut.genome, and cytobands will be the actual sections of introgression
  kpAddCytobandsAsLine(kp, color.table=c(acen = "#F92F27", stalk="#44DF94", gneg="#DADADA", gpos50="#626262", acen2="#A029FF"), lwd=12) # adds in the cytobands (which are our introgression), and colors them based on the colors we define
  kpAddBaseNumbers(kp, tick.dist = 10000000, add.units = T, minor.ticks = T, tick.len = 60, minor.tick.len = 20) # adds tick marks are scale for each chromosome. Size and spacing of each tick can be adjusted by the parameters here
  kpText(kp, data=peanut.genome[1], x=136000000, y=-1.5, labels = sample) # add label in the top right corner to specify in the plot which sample we are looking at
  kpRect(kp, y0=0, y1=-2, data=peanut.genome, border="black", lwd=1) # create a black border around each chromosome
  
  # copies the plot to save it as a jpeg file
  dev.copy(jpeg, filename=paste0("Introgression Plots/", sample, ".jpg"))
  dev.off()

  # show that the plotting process is fully complete
  print(paste(sample,"plotting is done!"))
}
