suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(data.table))
library(stringr)
library(optparse)

# TODO: convert all readlines calls to open the filepath as a file object first so that it doesn't read all teh data in at once
find_case_file <- function(case_files, network_id) {
  for (file in case_files) {
    if (str_detect(file, network_id)) {
      return(file)
    }
  }
}

generate_grid_colors <- function(network_grids, case_list) {
  color_vector <- c()
  for (grid in network_grids) {
    if (grid %in% case_list) {
      color_vector <- append(color_vector, "dark red")
    } else {
      color_vector <- append(color_vector, "dark gray")
    }
  }
  return(color_vector)
}

# Read the phenotype file if it is provided. We will return a
# hash env where the keys are ids and the values are the phenotype
# status (or variant they carry if that is the provided value)
load_pheno_file <- function(pheno_filepath) {
  # Lets first create the hash environment to key all of
  # our key value pairs
  pheno_hash <- new.env(hash = TRUE)

  # Lets make sure the pheno_filepath exists
  if (!file.exists(pheno_filepath)) {
    stop(paste0("ERROR: The file, ", pheno_filepath, ", does not exist. Terminating program..."), call. = FALSE)
  }

  con <- file(pheno_filepath, "r")
  on.exit(close(con))
  lines <- readLines(pheno_filepath, warn = FALSE)

  line_counter <- 1
  for (line in lines) {
    lines <- readLines(pheno_filepath, warn = FALSE)
    # This condition checks if we are at the end of the file
    if (length(line) == 0) {
      break
    }
    # We need to split each line
    split_line <- unlist(strsplit(trimws(line), "\t"))

    # We also need to make sure the line only has 2 elements like we expect
    if (length(split_line) != 2) {
      stop(paste0("ERROR: encountered a malformed line at line, ", line_counter, " : ", line))
    }

    # We also need to check and make sure that the individual
    # is not already in the hash. We are assuming each row in
    # the phenotype file is a unique id
    if (exists(split_line[1], envir = pheno_hash)) {
      stop("ERROR: Encountered a duplicate id in the phenotype file. Please ensure that each row in the phenotype file is a unique id")
    }

    pheno_hash[[split_line[1]]] <- split_line[2]
  }

  return(pheno_hash)
}

# We need to map the length column to a value that we can use to color the circos plot
format_ibd_data_length <- function(ibd_dataframe) {
  df <- ibd_dataframe %>%
    mutate(width = case_when(
      length < 3 ~ "grey",
      length < 10 ~ "blue",
      TRUE ~ "red"
    ))
  return(df)
}

define_id_colors <- function(ids, pheno_hash, case_list) {
  color_vector <- c()
  for (grid in network_grids) {
    if (grid %in% case_list) {
      color_vector <- append(color_vector, "dark red")
    } else {
      color_vector <- append(color_vector, "dark gray")
    }
  }
  return(color_vector)
}
# main function that will read through the DRIVE networks file to generate a circos plot for each network
process_network_file <- function(network_filepath, runtime_state) {
  if (!file.exists(network_filepath)) {
    stop(paste0("The network file, ", network_filepath, ", was not found on the system. Terminating program"))
  }
  # open the file for reading
  con <- file(network_filepath, "r")
  on.exit(close(con))

  for (line in lines) {
    lines <- readLines(network_filepath, warn = FALSE)

    # Exit the for loop if there is no value returned from readLines
    if (length(line) == 0) {
      break
    }
    split_line <- unlist(strsplit(trimws(line), "\t"))
    # we first need to handle the header

    # The we need to iterate through each network
    # If we provided a network id then we only need to process that and then break
    if (!is.null(runtime_state$network_id) && split_line[1] == runtime_state$network_id) {
      # add logic here
      break
    } else if (as.integer(split_line[2]) >= runtime_state$min_network_size) {
      # add logic to generate circos plot
    }
  }
}
##### Script technically will start running here #############
# Next lines until "parse_args()" deal with the CLI interface
option_list <- list(
  make_option(c("-n", "--network"),
    type = "character", default = NULL,
    help = "Input network file", metavar = "character"
  ),
  make_option(c("-i", "--ibd"),
    type = "character", default = NULL,
    help = "IBD file that list the pairwise IBD segments shared between individuals in the cohort. Currently the program only supports the hap-IBD format. This file should be tab separated and shouldn't have a header.", metavar = "character"
  ),
  make_option(c("-p", "--phenotype"),
    type = "character", default = NULL,
    help = "Tab separated text file that can be used to classify individuals with an affection status (This could be for a phenotype of interest or carriers of a variant. The file should have two columns: 'GRID' and 'Status', where each row is an id and then the corresponding status). If this file is provided then we will ignore the case id status in the DRIVE network file", metavar = "character"
  ),
  make_option("--id",
    type = "character", default = NULL,
    help = "Network ID from the DRIVE networks file. We will only generate plots for this network. If this id is not found in the file then the program will throw an error.", metavar = "character"
  ),
  make_option("--min-network-size",
    type = "integer", default = 2,
    help = "Minimum network size to plot (default %default)", metavar = "number"
  ),
  make_option("--pheno-column",
    type = "character", default = NULL,
    help = "Phenotype column name in the DRIVE network file (optional). This flag should not be provided if a phenotype file is provided.", metavar = "character"
  ),
  make_option(c("-o", "--output"),
    type = "character", default = "test.png",
    help = "Output image path for the circos plot. If the user is making a circos plot for every network in the file. Then the user should give a directory instead of a file name", metavar = "character"
  )
)

opt_parser <- OptionParser(
  option_list = option_list,
  description = "A tool to generate circos plots for visualizing IBD networks identified by DRIVE."
)
opt <- parse_args(opt_parser)

# TODO: Add someway to record what pararmeters are being
# used in the CLI

# This is a named list that we will use to keep track of things in
# our program
runtime_state <- list(
  network_id = opt$id,
  min_network_size = opt$`min-network-size`,
  pheno_col = opt$`pheno-column`
)

# If the user didn't pass any arguments then we need to print the help message
if (is.null(opt$network) || is.null(opt$ibd)) {
  print_help(opt_parser)
  stop("Both --network and --ibd arguments are required.\n", call. = FALSE)
}

# Check and make sure that the phenotype file and pheno-column arguments
# were not both provided
if (!is.null(opt$phenotype) && !is.null(opt$`pheno-column`)) {
  stop("ERROR: found values for both the '--phenotype' and '--pheno-column' flags. Program expects only 1 of these values to be provided")
}

# reading in the phenotype file if its provided
if (!is.null(opt$phenotype)) {
  print(paste0("Loading in the phenotype file: ", opt$phenotype))
  runtime_state$pheno_hash <- load_pheno_file(opt$phenotype)
} else {
  print(paste0("No phenotype file provided. Using the case status within the DRIVE file column: ", opt$`pheno-column`))
}

# Now we can read in the ibd_data. A dataframe will probably work
# best for this because we are going to access it by all of the ids
# in the network
print(paste0("Loading in the IBD segment data: ", opt$ibd))

runtime_state$ibd_df <- fread(opt$ibd, sep = "\t", header = FALSE)

colnames(runtime_state$ibd_df) <- c("pair_1", "hapID1", "pair_2", "hapID2", "chr", "start", "end", "length")

# Now lets filter our datatable for our cohort of interest.
# Only do this if the phenotype file was provided
if (runtime_state$pheno_hash != NULL) {
  cohort_ids <- ls(runtime_state$pheno_hash)

  runtime_state$ibd_df <- runtime_state$ibd_df[pair_1 %in% cohort_ids & pair_2 %in% cohort_ids]
}

## Now we can iterate through the DRIVE file

for (file in ibd_files) {
  split_filename <- unlist(strsplit(file, "_"))

  output_dir <- paste(split_filename[1], split_filename[2], "files", sep = "_")

  network_id <- paste(split_filename[4], split_filename[5], sep = "_")
  case_file <- find_case_file(case_files, network_id)

  ## Load the file into a dataframe

  ibd <- read.table(file, header = T, sep = "\t", stringsAsFactors = F)

  carriers <- read.table(case_file, stringsAsFactors = FALSE)

  colnames(carriers) <- c("grids")


  # ## setting the column names to certain values
  colnames(ibd) <- c("pair_1", "hapID1", "pair_2", "hapID2", "chr", "start", "end", "random column", "another random column", "length", "another column that usually equals 1")
  # print(colnames(ibd))
  color_list <- generate_grid_colors(unique(c(ibd$pair_1, ibd$pair_2)), carriers$grids)
  # ## creating a width column where if the value is greater than 10 then a 10 is entered in the column.
  # ## If the value is >3 but <10 then a 3 is entered in. If the value is <3 then a 1 is entered.
  ibd$width <- sapply(ibd$len, function(x) ifelse(x > 10, 10, 3))


  # ## if the value is greater than 10 then it gets assigned red. and if it is <10 but >3 then it is blue.
  # ## if it is <10 and <3 then it is colored gray
  # ibd$col = sapply(ibd$hapibd_len, function(x) ifelse(x>10, 'red', ifelse(x>3, 'blue', 'gray')))
  ibd$col <- sapply(ibd$length, function(x) ifelse(x > 10, "red", "blue"))

  # print(colnames(ibd))

  # ## This function will set the status as red for any individual identified as a carrier, and it will set the status to green if
  # ## the potential_missed_carrier was 1 and gray if it was 0
  # ibd$status = mapply(function(x, y) ifelse(x == 1, 'red', ifelse(y == 1, 'green', 'gray')), ibd$carrier_status, ibd$potential_missed_carrier)


  # ibd$status = sapply(function(x) ifelse(x %in% carriers$node, 'red', 'grey'), ibd$pair_1)

  # # This will basically create a list of red and grey or grey values of a length equivalent to the times each value shows up
  # # grid.col = c(rep('dark red', times=length(unique(ibd$pair_1))), rep('green', times=length(ibd[ibd$status == 'green',]$pair_2)), rep('dark gray', times=length(setdiff(ibd$pair_2, unique(ibd$pair_1)))))

  # grid.col = c(rep('dark red', times=length(carriers$grids)), rep('dark gray', times=length(setdiff(unique(c(ibd$pair_1,ibd$pair_2)), unique(carriers$grids)))))
  grid.col <- color_list

  # print(grid.col)
  # # setting the names of the list to the grids
  names(grid.col) <- c(unique(c(ibd$pair_1, ibd$pair_2)))
  # # names(grid.col) = c('R207504960','R208568059','R238563816','R251004118','R279769769', 'R256735849','R272528535','R280395674', 'R248254849','R201224925','R202872777','R201224925','R230165547', 'R223145375','R268938409', setdiff(unique(c(ibd$ID1, ibd$ID2)), c('R207504960','R208568059','R238563816','R251004118','R279769769', 'R256735849','R272528535','R280395674', 'R248254849','R201224925','R256735849','R223145375','R268938409','R202872777','R201224925','R230165547')))

  # grid.order = c(unique(carriers$grids), setdiff(unique(c(ibd$pair_1,ibd$pair_2)), unique(carriers$grids)))


  lgd_lines1 <- Legend(
    at = c("HCM Case", "HCM control"),
    legend_gp = gpar(fill = c("dark red", "dark gray")), labels_gp = gpar(fontsize = 15), title_position = "topleft",
    title = "Status"
  )

  lgd_lines2 <- Legend(
    at = c("> 10cM", "> 3cM"), type = "lines",
    legend_gp = gpar(col = c("red", "blue")), labels_gp = gpar(fontsize = 15), title_position = "topleft",
    title = "IBD Segment length"
  )

  lgd_list <- packLegend(lgd_lines1, lgd_lines2)

  output_name <- paste(output_dir, "/", network_id, "_circos_plot.png", sep = "")

  print(paste("writing output to", output_name, sep = " "))
  png(output_name, height = 12, width = 15, units = "in", res = 500)



  chordDiagram(ibd[, c("pair_1", "pair_2", "width")], col = ibd$col, grid.col = grid.col, annotationTrack = c("name", "grid"), annotationTrackHeight = c(0.02, 0.05), scale = T)

  # #order = grid.order
  # # highlight.sector(c('R207504960','R208568059','R238563816'),track.index = 1, col = "red",cex = 0.8, niceFacing = TRUE)
  draw(lgd_list, just = "left", x = unit(1.03, "snpc"), y = unit(0.15, "snpc"))

  dev.off()

  # print("finished creating plot")
}
