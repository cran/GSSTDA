## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(GSSTDA)

## -----------------------------------------------------------------------------
data("full_data")
data("survival_time")
data("survival_event")
data("case_tag")

## -----------------------------------------------------------------------------
# Gene selection information
gen_select_type <- "Top_Bot"
percent_gen_select <- 10 # Percentage of genes to be selected

## -----------------------------------------------------------------------------
#Mapper information
num_intervals <- 10
percent_overlap <- 40
distance_type <- "correlation"
clustering_type <- "hierarchical"
linkage_type <- "single" # only necessary if the type of clustering is hierarchical 
# num_bins_when_clustering <- 10 # only necessary if the type of clustering is hierarchical 
                                 # and the optimal_clustering_mode is "standard"
                                 # (this is not the case)

## -----------------------------------------------------------------------------
dsga_object <- dsga(full_data, survival_time, survival_event, case_tag)



## -----------------------------------------------------------------------------
gene_selection_object <- gene_selection(dsga_object, gen_select_type, percent_gen_select)


## -----------------------------------------------------------------------------
# Create data object
data_object <- list("full_data" = full_data, "survival_time" = survival_time,
                 "survival_event" = survival_event, "case_tag" = case_tag)
class(data_object) <- "data_object"


#Select gene from data object
gene_selection_object <- gene_selection(data_object, gen_select_type, percent_gen_select)



## -----------------------------------------------------------------------------
mapper_object <- mapper(data = gene_selection_object[["genes_disease_component"]], 
                        filter_values = gene_selection_object[["filter_values"]],
                        num_intervals = num_intervals,
                        percent_overlap = percent_overlap, distance_type = distance_type,
                        clustering_type = clustering_type,
                        linkage_type = linkage_type)



## -----------------------------------------------------------------------------
dsga_information <- results_dsga(dsga_object[["matrix_disease_component"]], case_tag)
print(dsga_information)

## -----------------------------------------------------------------------------
print(mapper_object)

## -----------------------------------------------------------------------------
plot_mapper(mapper_object)

## -----------------------------------------------------------------------------
gsstda_obj <- gsstda(full_data = full_data, survival_time = survival_time, 
                     survival_event = survival_event, case_tag = case_tag, 
                     gen_select_type = gen_select_type, 
                     percent_gen_select = percent_gen_select, 
                     num_intervals = num_intervals, 
                     percent_overlap = percent_overlap, 
                     distance_type = distance_type, 
                     clustering_type = clustering_type, 
                     linkage_type = linkage_type)




## -----------------------------------------------------------------------------
dsga_information <- results_dsga(gsstda_obj[["matrix_disease_component"]], case_tag)
print(dsga_information)

## -----------------------------------------------------------------------------
print(gsstda_obj[["mapper_obj"]])

## -----------------------------------------------------------------------------
plot_mapper(gsstda_obj[["mapper_obj"]])

