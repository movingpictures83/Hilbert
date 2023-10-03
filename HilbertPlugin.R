# ---------------------------------------------------------------------------------------------------------------------
#
#   hilbert.R (BETA)
#
#   This software is a "Camilo Valdes Work" under the terms of the United States Copyright Act.
#   Please cite the author(s) in any work or product based on this material. A prototype version of this script was
#   created as part of my doctoral work at the BioRG lab at Florida International University, and it was expanded
#   while at the Quantitative Life Sciences Initiative at the University of Nebraska-Lincoln.
#
#
#   OBJECTIVE:
#   The purpose of this program is to create Microbiome Maps using R. NOTE THAT THIS SOFTWARE IS A BETA VERSION.
#
#   AUTHOR:
#       Camilo Valdes
#
# ---------------------------------------------------------------------------------------------------------------------
library("IRanges")
library("circlize")
library("scales")
library("ComplexHeatmap")
library("RColorBrewer")
library("HilbertCurve")
require("optparse")

dyn.load(paste("RPluMA", .Platform$dynlib.ext, sep=""))
source("RPluMA.R")

options(width=512, digits=15, warn=1)

input <- function(inputfile) {
  parameters <<- read.table(inputfile, as.is=T);
  rownames(parameters) <<- parameters[,1];
  pfix = prefix()
  if (length(pfix) != 0) {
     pfix <- paste(pfix, "/", sep="")
  }

  
  ###############################################################
  # DEFAULT VALUES
  working_dir <<- paste(pfix, parameters["workdir", 2], sep="/")
  sample_id <<- parameters["sampleid", 2]
  col_name <<- 1
  col_abundance <<- 2
  verbose <<- FALSE
  labels2 <<- FALSE
  min_num_labels2 <<- 175  
  ###############################################################
  print(typeof(parameters))
  if ("colname" %in% rownames(parameters)) {
  col_name  <<- as.integer(parameters["colname", 2])
  }
  if ("colabundance" %in% rownames(parameters)) {
  col_abundance  <<- as.integer(parameters["colabundance", 2])
  }
  if ("verbose" %in% rownames(parameters)) {
 verbose  <<- as.boolean(parameters["verbose", 2])
  }
  if ("labels" %in% rownames(parameters)) {
  labels2   <<- as.boolean(parameters["labels", 2])
  }
  if ("minnumlabels" %in% rownames(parameters)) {
  min_num_labels2  <<- as.integer(parameters["minnumlabels", 2])
  }
  print("INPUT")
}

run <- function() {

}

output <- function(output_dir) {
print(paste( "[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] R Starting... ", sep=""))
    print(paste( "[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] Getting CLI Arguments...", sep=""))
    print(paste( "[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] Working DIR: ", working_dir, sep=""))
    print(paste( "[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] Output DIR: ", output_dir, sep=""))
    print(paste( "[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] Sample ID: ", sample_id, sep=""))
    print(paste( "[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] Column for Name: ", col_name, sep=""))
    print(paste( "[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] Column for Abundance: ", col_abundance, sep=""))
    print(paste( "[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] Image with Labels: ", labels2, sep=""))
    print(paste( "[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] Minimum No. of Children: ",
        min_num_labels2, sep=""))
    print(paste( "[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] Verbose Mode: ", verbose, sep=""))



working_dir = working_dir
print(working_dir)
#setwd(file.path( working_dir ))

base_analysis_dir = working_dir

display_labels2 = labels2
min_num_labels2 = min_num_labels2

col_microbe_name = col_name
col_abundance = col_abundance

output_dir = output_dir

sample_id = sample_id


# ------------------------------------------------------ Main ---------------------------------------------------------

print(paste("[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"]", sep=""))
print(paste("[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] Parsing Profile...", sep=""))

profiles_directory = base_analysis_dir

abundances_txt = paste(profiles_directory, "/", sample_id, ".txt", sep="")
abundances_table = read.delim(abundances_txt, header=FALSE, sep="\t", check.names=FALSE)

number_of_rows = nrow(abundances_table)
print(paste("[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] Number of Rows: ", number_of_rows, sep=""))

hilbert_level = 1
for(i in 1:20) {
    if(number_of_rows < 4^i){
       hilbert_level = i
       break
    }
}
print(paste("[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] Curve Level: ", hilbert_level, sep=""))

# ------------------------------------------------ Taxonomic Locales --------------------------------------------------

print(paste("[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"]", sep=""))
print(paste("[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] Parsing Taxonomic Locales...", sep=""))

previous_genus = ""
coord_prev = 1
coord_counter = 1
number_of_genera = 0

range_start_list = list()
range_end_list = list()

range_label_min_size = min_num_labels2
range_labels2_start_list = list()
range_labels2_end_list = list()
range_labels2_text = list()

abundances_table_sorted_by_genus = abundances_table[order(abundances_table[,col_microbe_name]),]
abundances_table_as_char = as.character(abundances_table_sorted_by_genus[,col_microbe_name])

for(bacteria_data in abundances_table_as_char) {
    bacteria_fields = strsplit(bacteria_data, '\t')
    bacteria_name = paste(unlist(bacteria_fields[1]), collapse=' ')
    strain_fields = strsplit(bacteria_name, ' ')
    genus = strain_fields[[1]][1]

    if(previous_genus == ""){
        previous_genus = genus
        number_of_genera = number_of_genera + 1
    }
    else if(previous_genus != genus){
        range_start = coord_prev
        range_end = coord_counter - 1
        range_size = range_end - range_start

        range_start_list = c(range_start_list, range_start)
        range_end_list = c(range_end_list, range_end)

        if(range_size >= range_label_min_size){
            range_labels2_start_list = c(range_labels2_start_list, range_start)
            range_labels2_end_list = c(range_labels2_end_list, range_end)
            range_labels2_text = c(range_labels2_text, previous_genus)
        }

        previous_genus = genus
        number_of_genera = number_of_genera + 1
        coord_prev = coord_counter

    } else {

    }
    coord_counter = coord_counter + 1
}

print(paste("[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] Unique Genera: ", number_of_genera, sep=""))
print(paste("[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"]", sep=""))

taxonomic_locales = IRanges(unlist(range_start_list), unlist(range_end_list))
taxonomic_locales_text = IRanges(unlist(range_labels2_start_list), unlist(range_labels2_end_list))

print(paste("[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] Starting Hilbert Image...", sep=""))

rescale_range_max = 1
rescale_range_min = 0

abundance_values    = log10(abundances_table_sorted_by_genus[,col_abundance])
abundance_values    = rescale(abundance_values, to=c(rescale_range_min, rescale_range_max))
abundance_max       = max(abundance_values[is.finite(abundance_values)])
abundance_median    = median(abundance_values[is.finite(abundance_values)])
abundance_min       = min(abundance_values[is.finite(abundance_values)])
number_of_taxa      = number_of_rows

print(paste("[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] Number of Taxa in Plot: ", number_of_rows, sep=""))
print(paste("[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] Abundance Max: ", abundance_max, sep=""))
print(paste("[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] Abundance Median: ", abundance_median, sep=""))
print(paste("[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] Abundance Min: ", abundance_min, sep=""))

abundance_ranges = c(abundance_min, abundance_median, abundance_max)

color_max = abundance_max
color_min = abundance_min
number_of_colors = 2

color_palette = c('white', 'red3')
color_function = colorRamp2(breaks=seq(color_min, color_max, length.out=number_of_colors), colors=color_palette)
legend_abundance = Legend(at=seq(color_min, color_max, by=number_of_colors/10), col_fun=color_function,
                          title="Abundance", direction="horizontal")
legend_list = (legend_abundance)

x = seq(1, number_of_rows, by=1)
x1 = x
x2 = x

# Draw the Hilbert Curve
scaling_factor = 8
figure_width_height = (2^hilbert_level) * scaling_factor

figure_png = paste(output_dir, "/", sample_id ,".png", sep="")
png(figure_png, width=figure_width_height, height=figure_width_height, res=72)

hilvert_curve = HilbertCurve(1, number_of_taxa, level=hilbert_level, mode="pixel", start_from="topleft")

#dev.off()

print(paste("[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"]", sep=""))
print(paste("[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"] Done.", sep=""))
print(paste("[", format(Sys.time(), "%m/%d/%y %H:%M:%S"),"]", sep=""))
}
