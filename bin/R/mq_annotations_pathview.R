#!/usr/bin/env Rscript
library("biomaRt")
library("optparse")
library(optparse)
library(yaml)
library(readr)
library(pathview)
library(xml2)
#install.packages("magick")
library(png)

# Define the command-line options
option_list <- list(
		      make_option(c("-c", "--config"), type = "character", default = NULL,
				                help = "Path to the config file", metavar = "character"),
		      
		      make_option(c("-o", "--outdir"), type = "character", default = ".",
				                help = "Output directory [default = %default]", metavar = "character")
		      )

# Parse the options
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Check if the config argument was provided
if (is.null(opt$config)) {
	  print_help(opt_parser)
  stop("Error: The config file must be specified using the -c or --config option.", call. = FALSE)
}

# Print out the options for debugging purposes
print(opt)

# Now you can use `opt$config` and `opt$outdir` in your script
config_path <- opt$config
output_directory <- opt$outdir

# Example: Print the paths
cat("Config file path:", config_path, "\n")
cat("Output directory:", output_directory, "\n")


#Read the config.yaml file
config <- yaml.load_file(config_path)

# Extract the output directory and reference list
references <- names(config$reference)

# List of colors in hexadecimal format
color_palette <- c(
	# Primary Colors
	"#FF0000",  # Red
	"#00FF00",  # Green
	"#0000FF",  # Blue

	# Secondary Colors
	"#FFFF00",  # Yellow
	"#FF00FF",  # Magenta
	"#00FFFF",  # Cyan

	# Additional Colors
	"#FFA500",  # Orange
	"#800080",  # Purple
	"#008000",  # Dark Green
	"#000080",  # Navy Blue
	"#FFC0CB",  # Pink
	"#A52A2A",  # Brown
	"#808080",  # Gray
	"#FFD700",  # Gold
	"#4B0082",  # Indigo
	"#00CED1",  # Dark Turquoise
	"#B0E0E6",  # Powder Blue
	"#2E8B57",  # Sea Green
	"#DC143C",  # Crimson
	"#4682B4",  # Steel Blue
	"#D2691E",  # Chocolate
	"#FF1493",  # Deep Pink
	"#9ACD32",  # Yellow Green
	"#40E0D0",  # Turquoise
	"#800000",  # Maroon
	"#8B4513",  # Saddle Brown
	"#5F9EA0",  # Cadet Blue
	"#8A2BE2",  # Blue Violet
	"#FF6347",  # Tomato
	"#FF4500",  # Orange Red
	"#7FFF00",  # Chartreuse
	"#DDA0DD",  # Plum
	"#98FB98",  # Pale Green
	"#F08080",  # Light Coral
	"#7B68EE",  # Medium Slate Blue
	"#FF69B4",  # Hot Pink
	"#E9967A",  # Dark Salmon
	"#BDB76B",  # Dark Khaki
	"#8FBC8F",  # Dark Sea Green
	"#556B2F",  # Dark Olive Green
	"#6495ED"   # Cornflower Blue
)



annotated_xml <- function(pathway_id, species, gene_data,  kegg.dir, legend=NULL, legend_names=NULL) {
	print(head(gene_data))
	
	ko_id <- gene_data[["Kegg Orthology ID"]]
	mapped_colors <- NULL
	mapping <- NULL
        label <- gene_data[['BLASTP']]
        #label <- paste(gene_data[['BLASTP']], "(", gene_data[['Gene Name']], ")", sep = "")
	k_id_to_protein <- setNames(label, ko_id)

	if (is.null(legend)){
		color <- color_palette[1]
	} else { 
		if (!is.null(legend_names)){
			sorted_values <- legend_names 
		} 
		else {
			sorted_values <- sort(unique(gene_data[[legend]]))
		}
		
		filt_values <- sort(unique(gene_data[[legend]]))
		print(head(gene_data))
		print(legend)
		print(gene_data[[legend]])
		# Map the sorted values to the color palette
		mapped_colors <- setNames(color_palette[1:length(sorted_values)], sorted_values)
		mapped_colors <- mapped_colors[filt_values]
		# Create a list to map each unique strain to its corresponding kid values
		# Create a named vector to map kid to strains
		print(mapped_colors)	
		legend_vals <- gene_data[[legend]]
		mapping <- setNames(legend_vals, ko_id)
		#mapping_list <- tapply(gene_data[["Kegg Orthology ID"]]  , gene_data[[legend]], function(x) x)
	}	
	xml_file_path= paste(kegg.dir, '/', 'ko', pathway_id, '.xml', sep='')
	xml_file_path_2= paste(kegg.dir, '/', species, pathway_id, '.xml', sep='')
	png_file_path= paste(kegg.dir, '/', species, pathway_id, '.png', sep='')
	png_file_path_copy= paste(kegg.dir, '/', species, pathway_id, '.copy.png', sep='')
	ko_id_unique <- unique(ko_id)

	print(xml_file_path)
	  if (!file.exists(png_file_path_copy)) {
		download.kegg(pathway.id = pathway_id, species = species, kegg.dir = kegg.dir,
		      file.type=c("xml", "png"))
		file_copied <- file.copy(from = png_file_path, to = png_file_path_copy, overwrite = TRUE)
		download.kegg(pathway.id = pathway_id, species = 'ko', kegg.dir = kegg.dir,
		      file.type=c("xml", "png"))

	  } else {
		message("xml already exists: ", xml_file_path)
	        file_copied <- file.copy(from = png_file_path_copy, to = png_file_path, overwrite = TRUE)
	  }
	#print(gene_data)
	doc <- read_xml(xml_file_path,as_html = TRUE)
	doc2 <- read_xml(xml_file_path_2,as_html = TRUE)
	img <- readPNG(png_file_path)
	img_width <- ncol(img)
	img_height <- nrow(img)
	# Set the resolution and dimensions
	output_width <- img_width  # or adjust as needed
	output_height <- img_height  # or adjust as needed
	output_res <- 300  # Set to 300 DPI or higher for better quality
	
	png(png_file_path, width = output_width, height = output_height, res = output_res)
	# Ensure the plotting area is large enough
	par(mar = c(0, 0, 0, 0), oma = c(0, 0, 0, 0))  # No margins and no outer margins

	# Display the image
	plot(1:2, type = "n", xlim = c(0, img_width), ylim = c(0, img_height), xlab = "", ylab = "", asp = 1, axes = FALSE, ann = FALSE, xaxs = "i", yaxs = "i")
	
	rasterImage(img, 0, 0, ncol(img), nrow(img))

	# Example: Extract all nodes (e.g., genes, compounds)
	nodes <- xml_find_all(doc, "//entry")
	# Extract node attributes (e.g., name, type, coordinates)
	color='red'
	for (node in nodes) {
		node_id <- xml_attr(node, "id")
	  	node_name <- xml_attr(node, "name")
		
		# Step 1: Extract KEGG IDs from the string
		kegg_ids <- unlist(strsplit(node_name, " "))

		# Step 2: Remove the "ko:" prefix from each ID
		kegg_ids_clean <- sub("^ko:", "", kegg_ids)
		#print(mapping)
		#print(kegg_ids_clean)
		if (!is.null(mapping)) { 
			color=NULL 
			legend_vals <- mapping[kegg_ids_clean]	
			# Loop through each element and print if not NA
			for (val in legend_vals) {
				if (!is.na(val)) {
					print(kegg_ids_clean)
					print(val)
					
					color <- mapped_colors[val]
					print(color)
    				}
			}
			#if (!is.na(legend_vals)) {
			#	print(legend_vals)
			#}
		}
		else {
			color='red'
		}
		# Step 3: Check if any of the cleaned KEGG IDs occur in the other list
		matching_ids <- kegg_ids_clean[kegg_ids_clean %in% ko_id_unique]
		# If you want to just check if there's any match (TRUE/FALSE)
		any_match <- any(kegg_ids_clean %in% ko_id_unique)
		if (!any_match == TRUE) { 
			next
			#color="#D3D3D3"
	        }
		lookup <- matching_ids[[1]]
		protein_name <- k_id_to_protein[ lookup]

		node_type <- xml_attr(node, "type")
	   	# Skip nodes of type "compound"
	        if (node_type == "compound") {
		        next
		 }	   
	    	# Extract coordinates (graphics element contains position and size)
	   	 graphics <- xml_find_first(node, "graphics")
	      	x <- as.numeric(xml_attr(graphics, "x"))
	      	y <- as.numeric(xml_attr(graphics, "y"))
	      	graphics_type <- xml_attr(graphics, "type")
	        if (graphics_type != "line"  ) {
		        next
		 }	   

	        width <- as.numeric(xml_attr(graphics, "width"))
	        height <- as.numeric(xml_attr(graphics, "height"))
	        fgcolor <- xml_attr(graphics, "fgcolor")
	        coords_str <- xml_attr(graphics, "coords")
		# Split the string into individual numbers
		coords <- as.numeric(unlist(strsplit(coords_str, ",")))
		
		# Convert to a matrix or data frame for easy access to pairs
		coords_matrix <- matrix(coords, ncol = 2, byrow = TRUE)

		# Print the coordinates as pairs
		#print(coords_matrix)
		
		# Adjust the y-coordinates because the origin (0,0) is at the top-left in images
		
		
		coords_matrix[, 2] <- img_height - coords_matrix[, 2]
		
		# Save the plot with annotations as a PNG
		#dev.off()		
		# Draw the points on the image
		lw=1.0
		#points(coords_matrix, col = "red", pch = 19, cex=lw/2)
		# Optionally, draw lines between the points
		# Extract the first point
		first_point <- coords_matrix[1, ]
		#print(first_point)
		#print(coords_matrix)
		# Extract the last point
		last_point <- coords_matrix[nrow(coords_matrix), ]
		#points(first_point, col = "red", pch = 19, cex=lw * 2)
		# Add the points
		find_middle_point <- function(points) {
	  		middle_index <- ceiling(nrow(points) / 2)
  			return(points[middle_index, ])
			} 
		middle_point <- find_middle_point(coords_matrix)

		for (i in 1:(nrow(coords_matrix) - 1)) {
			  segments(coords_matrix[i, 1], coords_matrix[i, 2], coords_matrix[i + 1, 1], coords_matrix[i + 1, 2], col = color, lwd = lw)
		}
		#color <- fgcolor
		#points(coords_matrix[1, 1], coords_matrix[1, 2], col = color, pch = 19, cex=lw/2)
		#points(coords_matrix[nrow(coords_matrix), 1], coords_matrix[nrow(coords_matrix), 2], col = color, pch = 19, cex=lw/2)
 	}	
	nodes <- xml_find_all(doc2, "//entry") # read the species specific compounds, to color them
	for (node in nodes) {
		node_id <- xml_attr(node, "id")
	  	node_name <- xml_attr(node, "name")
		node_type <- xml_attr(node, "type")
	   	# Skip nodes of type "compound"
		# Extract coordinates (graphics element contains position and size)
		 graphics <- xml_find_first(node, "graphics")
		x <- as.numeric(xml_attr(graphics, "x"))
		y <- as.numeric(xml_attr(graphics, "y"))
		y <- img_height - y
		graphics_type <- xml_attr(graphics, "type")
		if (node_type == 'compound') {
			#print(node_type)
			width <- as.numeric(xml_attr(graphics, "width"))
			height <- as.numeric(xml_attr(graphics, "height"))
			fgcolor <- xml_attr(graphics, "fgcolor")
			transparent_fg <- rgb(0, 0, 0, alpha = 0) 
			symbols(x, y, circles = width/2, inches = FALSE, add = TRUE, fg = fgcolor, bg=fgcolor)
			}
	}

	# Example: Add labels as the final step
	nodes <- xml_find_all(doc, "//entry")
	# Extract node attributes (e.g., name, type, coordinates)
	color='red'
	for (node in nodes) {
		node_id <- xml_attr(node, "id")
	  	node_name <- xml_attr(node, "name")
		
		# Step 1: Extract KEGG IDs from the string
		kegg_ids <- unlist(strsplit(node_name, " "))

		# Step 2: Remove the "ko:" prefix from each ID
		kegg_ids_clean <- sub("^ko:", "", kegg_ids)
		#print(mapping)
		#print(kegg_ids_clean)
		if (!is.null(mapping)) { 
			color=NULL 
			legend_vals <- mapping[kegg_ids_clean]	
			# Loop through each element and print if not NA
			for (val in legend_vals) {
				if (!is.na(val)) {
					print(kegg_ids_clean)
					print(val)
					
					color <- mapped_colors[val]
					print(color)
    				}
			}
		}
		else {
			color='red'
		}
		# Step 3: Check if any of the cleaned KEGG IDs occur in the other list
		matching_ids <- kegg_ids_clean[kegg_ids_clean %in% ko_id_unique]
		# If you want to just check if there's any match (TRUE/FALSE)
		any_match <- any(kegg_ids_clean %in% ko_id_unique)
		if (!any_match == TRUE) { 
			next
			#color="#D3D3D3"
	        }
		lookup <- matching_ids[[1]]
		protein_name <- k_id_to_protein[ lookup]

		node_type <- xml_attr(node, "type")
	   	# Skip nodes of type "compound"
	        if (node_type == "compound") {
		        next
		 }	   
	    	# Extract coordinates (graphics element contains position and size)
	   	graphics <- xml_find_first(node, "graphics")
	      	graphics_type <- xml_attr(graphics, "type")
	        if (graphics_type != "line"  ) {
		        next
		 }	   
	        coords_str <- xml_attr(graphics, "coords")
		# Split the string into individual numbers
		coords <- as.numeric(unlist(strsplit(coords_str, ",")))
		
		# Convert to a matrix or data frame for easy access to pairs
		coords_matrix <- matrix(coords, ncol = 2, byrow = TRUE)

		coords_matrix[, 2] <- img_height - coords_matrix[, 2]
		
		find_middle_point <- function(points) {
			dx <- diff(points[,1])
			dy <- diff(points[,2])
			distances <- sqrt(dx^2 + dy^2)
			cumulative_distances <- c(0, cumsum(distances))
			if (nrow(points) == 2) {
				    x_mid <- (points[1,1] + points[2,1]) / 2
				    y_mid <- (points[1,2] + points[2,2]) / 2
		                  } else {
			total_distance <- cumulative_distances[length(cumulative_distances)]
			mid_index <- which.min(abs(cumulative_distances - total_distance/2))

			x_mid <- points[mid_index,1]
			y_mid <- points[mid_index, 2]
		      }
	  		middle_index <- ceiling(nrow(points) / 2)
  			#return(points[middle_index, ])
  			return(c(x_mid, y_mid))
			} 
		middle_point <- find_middle_point(coords_matrix)
		# Draw the semi-transparent rectangle behind the text
		# Calculate the width and height of the text to create a background rectangle
		text_width <- strwidth(protein_name, cex=0.3)
		text_height <- strheight(protein_name, cex=0.3)
		yoffset=-10
		#rect(
		#       xleft = middle_point[1] - text_width / 2 - 0.01,
		#       ybottom = middle_point[2] + yoffset - text_height / 2 - 0.01,
		#       xright = middle_point[1] + text_width / 2 + 0.01,
		#       ytop = middle_point[2] + yoffset + text_height / 2 + 0.01,
		#       col = rgb(1, 1, 1, alpha=0.5), # white with 50% transparency
		#	     border = NA
		#		   )
		text(middle_point[1] , middle_point[2] + yoffset , labels=protein_name, cex=0.25, col=color)
	}
     	if (! is.null(mapping)) {
	legend("bottomright", legend =names(mapped_colors), col = mapped_colors, cex=0.5, pch = 15, title = "Strains", inset = c(0.05, 0.05))
		      }
	dev.off()		
}

split_and_expand <- function(df, column_to_split, separator = ";", new_index_name = NULL) {
	# Split the specified column by the separator into a list of vectors
	df[[column_to_split]] <- strsplit(as.character(df[[column_to_split]]), separator)

	# Expand the rows by unlisting the list
	expanded_df <- df[rep(1:nrow(df), sapply(df[[column_to_split]], length)), ]

	# Assign the unlisted elements to the expanded DataFrame
	expanded_df[[column_to_split]] <- unlist(df[[column_to_split]])

	# Optionally rename the row names (index)
	if (!is.null(new_index_name)) {
	rownames(expanded_df) <- seq_len(nrow(expanded_df))
	names(expanded_df)[which(names(expanded_df) == "row.names")] <- new_index_name
	}
      	data <- expanded_df
	data_filt <- data[!is.na(data$"KO Entrez"), ]
        data_filt <- data_filt[!duplicated(data_filt$"KO Entrez"), ]
	gene_data <- data_filt$"KO Entrez" 
        
	# Correctly formatted gene expression data using Entrez Gene IDs
	#gene_data <- c("885328" , "885367" , "885243" ,"885041" ,"885305" )  # INHA
  	#rownames(data_cleaned) <- gene_data
	#data_cleaned$val <- 1
	gene_vals <- rep(1, length(gene_data))
	gene_vals2 <- gene_data 
	# Remove any NA values in gene_data

	names(gene_vals2) <- gene_data 
	
	gene_data <- gene_vals2
	return(list(dataframe = data, list_element = gene_data))
}

# Iterate through each reference
for (reference in references) {
  print(reference)
  # Construct the path to the CSV file
  #csv_path <- file.path(output_directory, "proteingroups.ko.csv")
  csv_path <- file.path(output_directory, "annotations/export", reference, "protein2kegg.csv")
  out_path <- file.path(output_directory, "annotations/export", reference, 'pathview')
  dir.create(out_path, showWarnings = FALSE, recursive = TRUE)
  setwd(out_path)
  # Check if the file exists before trying to read it
  if (file.exists(csv_path)) {
	      # Read the CSV file
	      data <- read_csv(csv_path)
      
      # Print a message indicating which file is being processed
      cat("Processing file:", csv_path, "\n")
          
          # Display the first few rows of the data (for example)
          res <- split_and_expand(data, "Kegg Orthology ID", separator = ";")
	  data_cleaned <- res$dataframe
	  gene_vals <- res$list_element
	  csv_path <- file.path(output_directory, "annotations/export", reference, "protein2kegg.csv")
	  kegg.dir <- paste(output_directory, "/annotations/custom_kegg/", reference, sep='')
	  # Create the directory if it does not exist
	  if (!file.exists(kegg.dir)) {
		    dir.create(kegg.dir, recursive = TRUE)
	    message("Directory created: ", kegg.dir)
	  } else {
		    message("Directory already exists: ", kegg.dir)
	  }
	  if (length(gene_vals) == 0){
	    next
	  }
	  #print(gene_vals2)
	  
	  annotated_xml("01100", 'mtu', data_cleaned,  kegg.dir) 
	  pathview(gene.data = gene_vals, 
		   pathway.id = "01100",
		   species = "mtu",
		   plot.col.key=F,
		   discrete=list(gene=TRUE,
				 cpd=TRUE),
		   kegg.dir = kegg.dir,
		   out.dir= out_path,
		   out.suffix= 'protein_identifications',
		   gene.idtype="entrez",
		   same.layer=FALSE, 
		   map.null = FALSE,    # Ensures that only mapped nodes (genes) are displayed 
		   kegg.native = T,
	  low = "blue", mid = "grey", high = "red") 
	  
	  acetyl_path <- file.path(output_directory, "annotations/export", reference, "Combined_Proteins_Nterm_Acetylation.csv")
		
	   acetyl_data <- read_csv(acetyl_path)
	   acetyl_data <- merge(acetyl_data, data_cleaned, by='BLASTP')
           gene_data <- acetyl_data[!is.na(acetyl_data$"KO Entrez"), ]
           gene_data  <- gene_data [!duplicated(gene_data$"KO Entrez"), ]
	   gene_data <- gene_data$"KO Entrez" 
	   #gene_vals <- rep(1, length(gene_data))
	   gene_vals <- gene_data 
	   names(gene_vals) <- gene_data 
	   legend_names = acetyl_data[['Strains']] 
	   legend_names <- sort(unique(legend_names))
	   annotated_xml("01100", 'mtu', acetyl_data,  kegg.dir, legend="Strains", legend_names=legend_names) 

	   pathview(gene.data = gene_vals, 
		   pathway.id = "01100",
		   species = "mtu",
		   plot.col.key=F,
		   discrete=list(gene=TRUE,
				 cpd=TRUE),
		   kegg.dir = kegg.dir,
		   out.dir= out_path,
		   out.suffix= 'nterm_acetyl_all',
		   gene.idtype="entrez",
		   same.layer=FALSE, 
		   map.null = FALSE,    # Ensures that only mapped nodes (genes) are displayed 
		   kegg.native = T, low = "blue", mid = "grey", high = "red") 
	  acetyl_targets_path <- file.path(output_directory, "annotations/export", reference, "Combined_Proteins_Nterm_Acetylation_Targets.csv")
		
	   acetyl_data <- read_csv(acetyl_targets_path)
	   print(head(data_cleaned))
	   print(head(acetyl_data))
	   acetyl_data <- merge(acetyl_data, data_cleaned, by='BLASTP')
           print(head(acetyl_data))
	   gene_data <- acetyl_data[!is.na(acetyl_data$"KO Entrez"), ]
           gene_data  <- gene_data [!duplicated(gene_data$"KO Entrez"), ]
	   gene_data <- gene_data$"KO Entrez" 
	   #gene_vals <- rep(1, length(gene_data))
	   gene_vals <- gene_data 
	   names(gene_vals) <- gene_data 
	  
	   annotated_xml("01100", 'mtu', acetyl_data,  kegg.dir, legend="Strains", legend_names=legend_names) 

	   pathview(gene.data = gene_vals, 
		   pathway.id = "01100",
		   species = "mtu",
		   plot.col.key=F,
		   discrete=list(gene=TRUE,
				 cpd=TRUE),
		   kegg.dir = kegg.dir,
		   out.dir= out_path,
		   out.suffix= 'nterm_acetyl_targets',
		   gene.idtype="entrez",
		   same.layer=FALSE, 
		   map.null = FALSE,    # Ensures that only mapped nodes (genes) are displayed 
		   kegg.native = T, low = "blue", mid = "grey", high = "red") 

	  acetyl_diff_path <- file.path(output_directory, "annotations/export", reference, "Combined_Proteins_Nterm_Acetylation_Differences.csv")
		
	   acetyl_data <- read_csv(acetyl_diff_path)
	   acetyl_data <- merge(acetyl_data, data_cleaned, by='BLASTP')
           gene_data <- acetyl_data[!is.na(acetyl_data$"KO Entrez"), ]
           gene_data  <- gene_data [!duplicated(gene_data$"KO Entrez"), ]
	   gene_data <- gene_data$"KO Entrez" 
	   #gene_vals <- rep(1, length(gene_data))
	   gene_vals <- gene_data 
	   names(gene_vals) <- gene_data 
	  
	   annotated_xml("01100", 'mtu', acetyl_data,  kegg.dir, legend="Strains", legend_names=legend_names) 

	   pathview(gene.data = gene_vals, 
		   pathway.id = "01100",
		   species = "mtu",
		   plot.col.key=F,
		   discrete=list(gene=TRUE,
				 cpd=TRUE),
		   kegg.dir = kegg.dir,
		   out.dir= out_path,
		   out.suffix= 'nterm_acetyl_differences',
		   gene.idtype="entrez",
		   same.layer=FALSE, 
		   map.null = FALSE,    # Ensures that only mapped nodes (genes) are displayed 
		   kegg.native = T, low = "blue", mid = "grey", high = "red") 
	  
	   tss_diff_path <- file.path(output_directory, "annotations/export", reference, "Combined_Proteins_TSS_Differences.csv")
		
	   diff_data <- read_csv(tss_diff_path)
	   diff_data <- merge(diff_data, data_cleaned, by='BLASTP')
           gene_data <- diff_data[!is.na(diff_data$"KO Entrez"), ]
           gene_data  <- gene_data [!duplicated(gene_data$"KO Entrez"), ]
	   gene_data <- diff_data$"KO Entrez" 
	   #gene_vals <- rep(1, length(gene_data))
	   gene_vals <- gene_data 
	   names(gene_vals) <- gene_data 
	   legend_names = diff_data[['Comparisons']] 
	   legend_names <- sort(unique(legend_names))
	  
	   annotated_xml("01100", 'mtu', diff_data,  kegg.dir, legend="Comparisons", legend_names=legend_names) 

	   pathview(gene.data = gene_vals, 
		   pathway.id = "01100",
		   species = "mtu",
		   plot.col.key=F,
		   discrete=list(gene=TRUE,
				 cpd=TRUE),
		   kegg.dir = kegg.dir,
		   out.dir= out_path,
		   out.suffix= 'tss_differences',
		   gene.idtype="entrez",
		   same.layer=FALSE, 
		   map.null = FALSE,    # Ensures that only mapped nodes (genes) are displayed 
		   kegg.native = T, low = "blue", mid = "grey", high = "red") 
        } else {
		    cat("File not found:", csv_path, "\n")
	  }
}


