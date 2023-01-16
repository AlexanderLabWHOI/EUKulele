pacman::p_load(ggplot2,dplyr,data.table)

## FUNCTIONS ##

process_sub_eukulele <- function(eukulele_dir, label_out, all_results, salmon_dir,
                               salmon_format="num", is_protein=FALSE) {
    curr_dir = file.path(eukulele_dir,"taxonomy_estimation")
    for (curr_file in list.files(curr_dir)) {
        eukulele_estimates = data.frame(fread(file.path(curr_dir, curr_file), sep = "\t")) %>%
                dplyr::mutate(full_classification=case_when(!is.na(full_classification)~full_classification,
                                                            TRUE~"Unclassified")) %>%
                tidyr::separate(full_classification, sep = ";",
                                                     into = c("Domain","Supergroup","Phylum","Class","Order",
                                                              "Family","Genus","Species"))%>% 
                dplyr::mutate(Domain = unlist(lapply(as.character(EUKulele_out$Domain),trimws)))%>% 
                dplyr::mutate(Supergroup = unlist(lapply(as.character(EUKulele_out$Supergroup),trimws)))%>% 
                dplyr::mutate(Phylum = unlist(lapply(as.character(EUKulele_out$Phylum),trimws)))%>% 
                dplyr::mutate(Class = unlist(lapply(as.character(EUKulele_out$Class),trimws)))%>% 
                dplyr::mutate(Order = unlist(lapply(as.character(EUKulele_out$Order),trimws)))%>% 
                dplyr::mutate(Family = unlist(lapply(as.character(EUKulele_out$Family),trimws)))%>% 
                dplyr::mutate(Genus = unlist(lapply(as.character(EUKulele_out$Genus),trimws)))%>% 
                dplyr::mutate(Species = unlist(lapply(as.character(EUKulele_out$Species),trimws)))
        
        if (is_protein) {
            eukulele_estimates = eukulele_estimates %>% tidyr::separate(transcript_name,sep="\\.p",
                                                                        into=c("transcript_name",
                                                                               "protein_id"))
        }
        number_file = unlist(strsplit(curr_file,"_"))[1]
        salmon_file = read.csv(file.path(salmon_dir,paste0(salmon_format,as.character(number_file), "_quant"),
                                 "quant.sf"), sep = "\t")
        matched_file= salmon_file %>% dplyr::full_join(eukulele_estimates,by=c("Name"="transcript_name"))
        matched_file["Sample"] = curr_file
        matched_file["SplitSamp"] = number_file
        matched_file["Type"] = label_out

        matched_file =  matched_file %>%
            tidyr::replace_na(list("Domain"="Uncertain",
                                   "Supergroup"="Uncertain",
                                   "Phylum"="Uncertain",
                                   "Class"="Uncertain",
                                   "Order"="Uncertain",
                                   "Family"="Uncertain",
                                   "Genus"="Uncertain",
                                   "Species"="Uncertain"))
        all_results = all_results %>% dplyr::bind_rows(matched_file)
    }
    return(all_results)
}

all_results=data.frame()
all_results = process_sub_eukulele("CAG_eukulele_pleuro_metazoans", "Metazoan_DB", all_results, 
                     "/vortexfs1/omics/alexander/akrinos/2021-remodeling-eukrhythmic/2021-09-ALOHA/intermediate-files/04-compare/09-CAG-mapping/salmon",
                     salmon_format="num", is_protein=TRUE)
