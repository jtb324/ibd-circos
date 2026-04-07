suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
library(stringr)
library(optparse)

#' Load Phenotype File
#'
#' Reads a tab-separated phenotype file and returns an environment (hash map)
#' where keys are individual IDs and values are their phenotype status.
#'
#' @param pheno_filepath Character string. Path to the phenotype file.
#'
#' @return An environment containing ID-phenotype key-value pairs.
#' @export
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

  line_counter <- 1
  while (TRUE) {
    lines <- readLines(con, n = 1, warn = FALSE)
    # This condition checks if we are at the end of the file
    if (length(lines) == 0) {
      break
    }
    # We need to split each line
    split_line <- unlist(strsplit(trimws(lines), "\t"))

    # We also need to make sure the line only has 2 elements like we expect
    if (length(split_line) != 2) {
      stop(paste0("ERROR: encountered a malformed line at line, ", line_counter, " : ", lines))
    }

    # We also need to check and make sure that the individual
    # is not already in the hash. We are assuming each row in
    # the phenotype file is a unique id
    if (exists(split_line[1], envir = pheno_hash)) {
      stop("ERROR: Encountered a duplicate id in the phenotype file. Please ensure that each row in the phenotype file is a unique id")
    }

    pheno_hash[[split_line[1]]] <- split_line[2]
    line_counter <- line_counter + 1
  }

  return(pheno_hash)
}

# Reading in the mapping file so that we can map ids in the DRIVE file. 
# This is only needed if our ids in the driv file are different the the 
# ids used in the project. We will return the hash were the keys are the 
# original id and the values are the mapped ids
load_map_file <- function(mapping_file) {
  map_hash <- new.env(hash = TRUE)
  reverse_hash_map <- new.env(hash = TRUE)

  if (!file.exists(mapping_file)) {
    stop(paste0("ERROR: The provided mapping file, ", mapping_file, ", does not exist. Terminating program..."), call. = FALSE)
  }
  # On the file for streaming and then set on.exit so that the file closes at the end of the function
  con <- file(mapping_file, "r")
  on.exit(close(con))

  line_counter <- 1
  # This is basically a do-while loop
  repeat {
    lines <- readLines(con, n = 1, warn = FALSE)

    split_line <- unlist(strsplit(trimws(lines), "\t"))

    # if the line is empty then we have hit the end of the file 
    # and need to break out of the loop
    if (length(lines) == 0) {
      break
    }

    if (length(split_line) != 2) {
      stop(paste0("ERROR: encountered a malformed line at line, ", line_counter, " : ", lines, ". Expected the first column to be the original ID and the second colum to be the ids that we wish to map values to"))
    }

    map_hash[[split_line[1]]] <- split_line[2]
    reverse_hash_map[[split_line[2]]] <- split_line[1]

    line_counter <- line_counter + 1

  }
  return(list(forward_map=map_hash, reverse_map=reverse_hash_map))
}

# We need to map the length column to a value that we can use to color the circos plot
#' Format IBD Data Length
#'
#' Maps IBD segment lengths to specific colors and widths for visualization.
#' Segments < 3cM are grey, < 10cM are blue, and >= 10cM are red.
#'
#' @param ibd_dataframe A data frame containing IBD segment information, 
#'   including a `length` column.
#'
#' @return A data frame with additional `col` and `width` columns.
#' @export
format_ibd_data_length <- function(ibd_dataframe) {
  df <- ibd_dataframe %>%
    mutate(
      col = case_when(
        length < 3 ~ "grey",
        length < 10 ~ "blue",
        TRUE ~ "red"
      ),
      width = length
    )
  return(df)
}

#' Get Phenotype Column Index
#'
#' Finds the index of a specific phenotype column name within a header line.
#'
#' @param header_line Character vector. The split header line of a file.
#' @param pheno_col_name Character string. The name of the phenotype column to find.
#'
#' @return Integer. The index of the phenotype column.
#' @export
get_case_col_indx <- function(header_line, pheno_col_name) {
  if (length(header_line) <= 1) {
    stop("ERROR: Encountered a malformed header line. Terminating program...", call. = FALSE)
  }
  col_indx <- which(header_line == pheno_col_name)

  if (length(col_indx) == 0) {
    stop(paste0("ERROR: unable to find the column, ", pheno_col_name, ", in the DRIVE results file. This probably indicates a typo. Terminating program...", call. = FALSE))
  }
  return(col_indx)
}

#' Generate Grid Colors
#'
#' Assigns colors to network members (grids) based on their phenotype status
#' from a hash or a provided list of cases.
#'
#' @param network_grids Character vector. IDs of individuals in the network.
#' @param pheno_hash Optional environment. Hash map of IDs to phenotypes.
#' @param case_list Optional character vector. List of IDs considered cases.
#'
#' @return A named character vector of colors.
#' @export
generate_grid_colors <- function(network_grids, pheno_hash = NULL, case_list = NULL) {
  color_vector <- sapply(network_grids, function(id) {
    if (!is.null(pheno_hash)) {
      status <- pheno_hash[[id]]
      if (is.null(status)) {
        return("dark gray")
      }
      return("dark red")
    } else if (!is.null(case_list)) {
      return(ifelse(id %in% case_list, "dark red", "dark gray"))
    } else {
      return("dark gray")
    }
  })
  names(color_vector) <- network_grids
  return(color_vector)
}

#' Generate Circos Plots
#'
#' Creates and saves a circos plot for a specific network, showing IBD connections
#' between members.
#'
#' @param network_id Character string. The identifier for the network.
#' @param network_members Character vector. IDs of individuals in the network.
#' @param cases Character vector. IDs of individuals considered cases in this network.
#' @param runtime_state List. Global state containing IBD data and output configuration.
#'
#' @return None. Saves a PNG file to the output directory.
#' @export
generate_circos_plots <- function(network_id, network_members, cases, runtime_state) {
  # Filter IBD data for pairs where BOTH are in the network
  network_segments <- runtime_state$ibd_df[pair_1 %in% network_members & pair_2 %in% network_members]

  if (!is.null(runtime_state$mapped_id_list)) {
    network_segments$pair_1 <- sapply(network_segments$pair_1, function(x) get(x, envir=runtime_state$map_hash))
    network_segments$pair_2 <- sapply(network_segments$pair_2, function(x) get(x, envir=runtime_state$map_hash))
  }

  if (nrow(network_segments) == 0) {
    print(paste("No IBD segments found for network", network_id))
    return()
  }

  if (runtime_state$save_segments) {
    segments_output_name <- paste(runtime_state$output, "/", network_id, "_ibd_segments.txt", sep = "")
    print(paste("Saving segments to", segments_output_name, sep = " "))
    fwrite(network_segments, segments_output_name, sep = "\t")
  }

  # We need to map the ids if we have provided a hash
  if (!is.null(runtime_state$map_hash)) {
    mapped_ids = map_ids(network_members, runtime_state$map_hash)
  } else {
    mapped_ids = network_id
  }

  grid.col <- generate_grid_colors(mapped_ids, runtime_state$pheno_hash, cases)

  # Add a column "col" and "width" that tell what color the chords
  # are in the chart
  network_segments <- format_ibd_data_length(network_segments)

  # If the user provides either a phenotype file or determines cases from the drive file 
  # then we need to add a upper legend that has the case/control colors. Otherwise, we only 
  # need to generate the legend for the colors based on IBD segment length. The following 
  # lines construct a list of legend objects and then passes that to the pack legend call 
  # that will dynamically construct the legend
  lgd_obj_list = list()

  if (!is.null(runtime_state$pheno_hash) || !is.null(cases)) {
    lgd_lines1 <- Legend(
      at = c("Case", "control"),
      legend_gp = gpar(fill = c("dark red", "dark gray")), labels_gp = gpar(fontsize = 15), title_position = "topleft",
      title = "Status"
    )
    lgd_obj_list[[length(lgd_obj_list)+1]] <- lgd_lines1
  } 

  lgd_lines2 <- Legend(
    at = c("> 10cM", "> 3cM", "> 1cM"), type = "lines",
    legend_gp = gpar(col = c("red", "blue", "grey")), labels_gp = gpar(fontsize = 15), title_position = "topleft",
    title = "IBD Segment length"
  )

  lgd_obj_list[[length(lgd_obj_list)+1]] <- lgd_lines2

  lgd_list <- do.call(packLegend, lgd_obj_list)

  output_name <- paste(runtime_state$output, "/", network_id, "_circos_plot.png", sep = "")

  print(paste("writing output to", output_name, sep = " "))
  png(output_name, height = 12, width = 15, units = "in", res = 500)

  # chordDiagram expects a dataframe with [from, to, value] and optional colors
  circos.clear()
  chordDiagram(network_segments[, c("pair_1", "pair_2", "width")],
    col = network_segments$col,
    grid.col = grid.col,
    annotationTrack = c("name", "grid"),
    annotationTrackHeight = c(0.02, 0.05),
    scale = TRUE
  )

  draw(lgd_list, just = "left", x = unit(1.03, "snpc"), y = unit(0.15, "snpc"))

  dev.off()
}

# Map our ids that we have read in to the new ids. Function takes an array of the current ids and then maps 
# each id to a new value and returns that array. We return the new array
map_ids <- function(current_ids, mapping_hash) {
  new_array <- sapply(current_ids, function(x, h = mapping_hash) {
    if (exists(x, envir = h)) {
      return(h[[x]])
    } else {
      stop(paste0("Missing a mapping id for participant, ", x, ". Please make sure every id has a mapping"))
    }
  })
  return(new_array)
}

# main function that will read through the DRIVE networks file to generate a circos plot for each network
#' Process Network File
#'
#' Iterates through a DRIVE network file and generates circos plots for qualifying networks.
#'
#' @param network_filepath Character string. Path to the DRIVE network results file.
#' @param runtime_state List. Global state containing IBD data and configuration.
#'
#' @return None.
#' @export
process_network_file <- function(network_filepath, runtime_state) {
  if (!file.exists(network_filepath)) {
    stop(paste0("The network file, ", network_filepath, ", was not found on the system. Terminating program"))
  }
  # open the file for reading
  con <- file(network_filepath, "r")
  on.exit(close(con))

  # this variable is how we keep track of whether we have processed the
  # header line.
  header_processed <- FALSE
  pheno_col_indx <- NULL

  while (TRUE) {
    lines <- readLines(con, n = 1, warn = FALSE)

    # Exit the loop if there is no value returned from readLines
    if (length(lines) == 0) {
      break
    }
    split_line <- unlist(strsplit(trimws(lines), "\t"))

    # we first need to handle the header
    if (!header_processed) {
      if (is.null(runtime_state$pheno_hash) && !is.null(runtime_state$pheno_col)) {
        pheno_col_indx <- get_case_col_indx(split_line, runtime_state$pheno_col)
      }
      header_processed <- TRUE
      next
    }

    if (!is.null(pheno_col_indx)) {
      cases <- unlist(strsplit(split_line[pheno_col_indx], ","))
    } else {
      cases <- NULL
    }

    network_id <- split_line[1]
    network_size <- as.integer(split_line[2])
    network_members <- unlist(strsplit(split_line[7], ","))

    # If we have loaded in the mappings then we need to update the network 
    # members list here. The order of the mapped ids will be the same as 
    # the original list
    if (!is.null(runtime_state$map_hash)) {
      mapped_id_list = map_ids(network_members, runtime_state$map_hash)
      runtime_state$mapped_id_list = mapped_id_list
    } else {
      runtime_state$mapped_id_list = NULL
    }

    # If we provided a network id then we only need to process that and 
    # then break
    if (!is.null(runtime_state$network_id)) {
      if (network_id == runtime_state$network_id) {
        generate_circos_plots(network_id, network_members, cases, runtime_state)
        break
      }
    } else {
      # We also need to make sure that the networks are above the size threshold
      if (network_size >= runtime_state$min_network_size) {
        generate_circos_plots(network_id, network_members, cases, runtime_state)
      }
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
    help = "Tab separated text file that can be used to classify individuals with an affection status. The file should have two columns: 'GRID' and 'Status'.", metavar = "character"
  ),
  make_option(c("-m", "--map-file"),
    type = "character", default = NULL,
    help = "Filepath to a tab separated text file that maps our ids to new values. We expect 2 columns where the first column is the original id and the second column is the new id. We also expect a header", metavar = "character"
  ),
  make_option("--id",
    type = "character", default = NULL,
    help = "Network ID from the DRIVE networks file.", metavar = "character"
  ),
  make_option("--min-network-size",
    type = "integer", default = 2,
    help = "Minimum network size to plot (default %default)", metavar = "number"
  ),
  make_option("--pheno-column",
    type = "character", default = NULL,
    help = "Phenotype column name in the DRIVE network file (optional).", metavar = "character"
  ),
  make_option("--save-segments",
    action = "store_true", default = FALSE,
    help = "Save the IBD segments used for plotting each individual circos plot (default %default)", metavar = "boolean"
  ),
  make_option(c("-o", "--output"),
    type = "character", default = "test/",
    help = "Output directory that the circos plot will be written to. (default %default)", metavar = "character"
  )
)

opt_parser <- OptionParser(
  option_list = option_list,
  description = "A tool to generate circos plots for visualizing IBD networks identified by DRIVE."
)
opt <- parse_args(opt_parser)

# TODO: Add someway to record what pararmeters are being
# used in the CLI
# This is a named list that we will use to keep track of things
runtime_state <- list(
  network_id = opt$id,
  min_network_size = opt$`min-network-size`,
  pheno_col = opt$`pheno-column`,
  save_segments = opt$`save-segments`,
  output = opt$output
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

if (!dir.exists(opt$output)) {
  print(paste0("The output directory, ", opt$output, ", was not found. Creating the output directory"))
  dir.create(opt$output, recursive = TRUE)
}

# reading in the phenotype file if its provided
if (!is.null(opt$phenotype)) {
  print(paste0("Loading in the phenotype file: ", opt$phenotype))
  runtime_state$pheno_hash <- load_pheno_file(opt$phenotype)
} else {
  print(paste0("No phenotype file provided. Using the case status within the DRIVE file column: ", opt$`pheno-column`))
  runtime_state$pheno_hash <- NULL
}

# If we have the mapping file provided then we need to read that in and 
# store it in our state. Otherwise we will just set the attribute to NULL
if (!is.null(opt$`map-file`)) {
  print(paste0("Loading in mappings from the mapping file: ", opt$`map-file`))
  return_val <- load_map_file(opt$`map-file`)

  runtime_state$map_hash = return_val$forward_map
  runtime_state$reverse_map_hash = return_val$reverse_map
} else {
  runtime_state$map_hash <- NULL 
  runtime_state$reverse_map_hash <- NULL
}

# Now we can read in the ibd_data.
print(paste0("Loading in the IBD segment data: ", opt$ibd))

# fread can handle .gz files
# runtime_state$ibd_df <- fread(opt$ibd, sep = "\t", header = FALSE)
runtime_state$ibd_df <- fread(opt$ibd, sep = "\t", 
            header = FALSE, 
            col.names = c("pair_1", "hapID1", "pair_2", "hapID2", "chr", "start", "end", "length"),
            colClasses=list(character = c(1,3)))

# colnames(runtime_state$ibd_df) <- c("pair_1", "hapID1", "pair_2", "hapID2", "chr", "start", "end", "length")

print("Filtering the ibd segments to our cohort of interest")
# Now lets filter our datatable for our cohort of interest.
# Only do this if the phenotype file was provided
if (!is.null(runtime_state$pheno_hash)) {
  cohort_ids <- ls(runtime_state$pheno_hash)
  if (!is.null(runtime_state$map_hash)) {
    cohort_ids = map_ids(cohort_ids, runtime_state$reverse_map_hash)
  }
  runtime_state$ibd_df <- runtime_state$ibd_df[pair_1 %in% cohort_ids & pair_2 %in% cohort_ids]
}

## Now we can iterate through the DRIVE file
process_network_file(opt$network, runtime_state)

print("finished creating plots")
