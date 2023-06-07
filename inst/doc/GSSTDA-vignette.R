## ---- include = FALSE---------------------------------------------------------
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
num_intervals <- 5
percent_overlap <- 40
distance_type <- "cor"
clustering_type <- "hierarchical"
linkage_type <- "single" # only necessary if the type of clustering is hierarchical 
# num_bins_when_clustering <- 10 # only necessary if the type of clustering is hierarchical 
                                 # and the optimal_clustering_mode is "silhouette"
                                 # (this is not the case)

## -----------------------------------------------------------------------------
DGSA_object <- DGSA(full_data, survival_time, survival_event, case_tag)



## -----------------------------------------------------------------------------
geneSelection_object <- geneSelection(DGSA_object, gen_select_type, percent_gen_select)


## -----------------------------------------------------------------------------
# Create data object
data_object <- list("full_data" = full_data, "survival_time" = survival_time,
                 "survival_event" = survival_event, "case_tag" = case_tag)
class(data_object) <- "data_object"


#Select gene from data object
geneSelection_object <- geneSelection(data_object, gen_select_type, percent_gen_select)



## -----------------------------------------------------------------------------
mapper_object <- mapper(full_data = geneSelection_object[["genes_disease_component"]], 
                        filter_values = geneSelection_object[["filter_values"]],
                        num_intervals = num_intervals,
                        percent_overlap = percent_overlap, distance_type = distance_type,
                        clustering_type = clustering_type,
                        linkage_type = linkage_type, 
                        optimal_clustering_mode = optimal_clustering_mode)



## -----------------------------------------------------------------------------
DGSA_information <- results_DGSA(DGSA_object[["matrix_disease_component"]], case_tag)
print(DGSA_information)

## -----------------------------------------------------------------------------
print(mapper_object)

## -----------------------------------------------------------------------------
plot_mapper(mapper_object)

## -----------------------------------------------------------------------------
GSSTDA_obj <- GSSTDA(full_data = full_data, survival_time = survival_time, 
                     survival_event = survival_event, case_tag = case_tag, 
                     gen_select_type = gen_select_type, 
                     percent_gen_select = percent_gen_select, 
                     num_intervals = num_intervals, 
                     percent_overlap = percent_overlap, 
                     distance_type = distance_type, 
                     clustering_type = clustering_type, 
                     linkage_type = linkage_type)




## -----------------------------------------------------------------------------
DGSA_information <- results_DGSA(GSSTDA_obj[["matrix_disease_component"]], case_tag)
print(DGSA_information)

## -----------------------------------------------------------------------------
print(GSSTDA_obj[["mapper_obj"]])

## -----------------------------------------------------------------------------
plot_mapper(GSSTDA_obj[["mapper_obj"]])

