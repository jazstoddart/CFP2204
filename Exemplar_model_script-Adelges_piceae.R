require(dplyr)
 
    require(ape)

    require(brms)

    require(tidyr)

    require(ggplot2)

    require(distributional)

    remotes::install_github("stan-dev/cmdstanr")

    require(cmdstanr)

    install_cmdstan(cores=8,overwrite = TRUE)

    set_cmdstan_path("/data/home/mpx605/.cmdstan/cmdstan-2.33.1")

    

    ## Import data

    # Phloygenetic tree

    megatree <-read.tree("./Data/accepted_minitree_EURO_FI.tre")

    # Host status matrix (HSM) dataframe

    host_matrix_full<-read.csv("./Data/CABI_HostMatrix.csv")

    host_matrix_NA_full<-read.csv("./Data/CABI_HostMatrix_NA.csv")

    

    ## Specify the cluster of interest

    cluster<-c("Adelges_piceae",
"Megastigmus_rafni",
"Megastigmus_suspectus",
"Fomes_fomentarius",
"Cydia_splendana",
"Saperda_scalaris",
"Caliroa_annulipes",
"Tortrix_viridana",
"Xiphydria_camelus",
"Scolytus_scolytus",
"Paranthrene_tabaniformis",
"Fomitopsis_pinicola",
"Xyleborus_dispar",
"Heterobasidion_annosum",
"Saturnia_pyri",
"Neonectria_ditissima",
"Armillaria_ostoyae",
"Hylurgops_palliatus",
"Chionaspis_salicis",
"Heterobasidion_parviporum",
"Pityogenes_chalcographus",
"Sirococcus_conigenus",
"Gremmeniella_abietina",
"Hylastes_ater",
"Mycosphaerella_pini",
"Sphaeropsis_sapinea",
"Pissodes_castaneus")
## Cutting down the HSM to only include rows for plants with host status 

    host_matrix_temp<-cbind(host_matrix_full[,1],host_matrix_full[, colnames(host_matrix_full) %in% cluster ==TRUE])

    names(host_matrix_temp)[1] <- "Host"

    host_matrix<-host_matrix_temp[rowSums(host_matrix_temp[, -1]) > 0, ]

    host_matrix<-rbind(host_matrix,host_matrix_temp[(length(host_matrix_temp$Host)-27):length(host_matrix_temp$Host),])

    host_matrix_temp<-cbind(host_matrix_NA_full[,1],host_matrix_NA_full[, colnames(host_matrix_NA_full) %in% cluster ==TRUE])

    names(host_matrix_temp)[1] <- "Host"

    host_matrix_NA<-host_matrix

    host_matrix_NA[host_matrix_NA == 0] <- NA

    

    ## Generate phylogenetic covariance matrix

    A<-ape::vcv.phylo(megatree)

    A[1:length(megatree$tip.label),1:length(megatree$tip.label)]<-

      sqrt(A[1:length(megatree$tip.label),1:length(megatree$tip.label)])

    ## Set model form, cluster names without marks

    fullform <- bf(mvbind(
Adelges_piceae,
Megastigmus_rafni,
Megastigmus_suspectus,
Fomes_fomentarius,
Cydia_splendana,
Saperda_scalaris,
Caliroa_annulipes,
Tortrix_viridana,
Xiphydria_camelus,
Scolytus_scolytus,
Paranthrene_tabaniformis,
Fomitopsis_pinicola,
Xyleborus_dispar,
Heterobasidion_annosum,
Saturnia_pyri,
Neonectria_ditissima,
Armillaria_ostoyae,
Hylurgops_palliatus,
Chionaspis_salicis,
Heterobasidion_parviporum,
Pityogenes_chalcographus,
Sirococcus_conigenus,
Gremmeniella_abietina,
Hylastes_ater,
Mycosphaerella_pini,
Sphaeropsis_sapinea,
Pissodes_castaneus)
 ~ 1 + (1|p|gr(Host, cov = A)),
    family = bernoulli())
priors <- c(set_prior("student_t(2,5.51,1)", class = "Intercept", resp = "Adelgespiceae"),
set_prior("student_t(2,1.77,1)", class = "Intercept", resp = "Megastigmusrafni"),
set_prior("student_t(2,4.94,1)", class = "Intercept", resp = "Megastigmussuspectus"),
set_prior("student_t(2,3.94,1)", class = "Intercept", resp = "Fomesfomentarius"),
set_prior("student_t(2,4.53,1)", class = "Intercept", resp = "Cydiasplendana"),
set_prior("student_t(2,3.71,1)", class = "Intercept", resp = "Saperdascalaris"),
set_prior("student_t(2,1.24,1)", class = "Intercept", resp = "Caliroaannulipes"),
set_prior("student_t(2,3.72,1)", class = "Intercept", resp = "Tortrixviridana"),
set_prior("student_t(2,0.79,1)", class = "Intercept", resp = "Xiphydriacamelus"),
set_prior("student_t(2,2.64,1)", class = "Intercept", resp = "Scolytusscolytus"),
set_prior("student_t(2,4.19,1)", class = "Intercept", resp = "Paranthrenetabaniformis"),
set_prior("student_t(2,3.03,1)", class = "Intercept", resp = "Fomitopsispinicola"),
set_prior("student_t(2,5.5,1)", class = "Intercept", resp = "Xyleborusdispar"),
set_prior("student_t(2,4.93,1)", class = "Intercept", resp = "Heterobasidionannosum"),
set_prior("student_t(2,2.89,1)", class = "Intercept", resp = "Saturniapyri"),
set_prior("student_t(2,3.5,1)", class = "Intercept", resp = "Neonectriaditissima"),
set_prior("student_t(2,3.71,1)", class = "Intercept", resp = "Armillariaostoyae"),
set_prior("student_t(2,3.95,1)", class = "Intercept", resp = "Hylurgopspalliatus"),
set_prior("student_t(2,3.33,1)", class = "Intercept", resp = "Chionaspissalicis"),
set_prior("student_t(2,5.51,1)", class = "Intercept", resp = "Heterobasidionparviporum"),
set_prior("student_t(2,3.03,1)", class = "Intercept", resp = "Pityogeneschalcographus"),
set_prior("student_t(2,3.94,1)", class = "Intercept", resp = "Sirococcusconigenus"),
set_prior("student_t(2,4.2,1)", class = "Intercept", resp = "Gremmeniellaabietina"),
set_prior("student_t(2,2.89,1)", class = "Intercept", resp = "Hylastesater"),
set_prior("student_t(2,3.17,1)", class = "Intercept", resp = "Mycosphaerellapini"),
set_prior("student_t(2,4.52,1)", class = "Intercept", resp = "Sphaeropsissapinea"),
set_prior("student_t(2,1.43,1)", class = "Intercept", resp = "Pissodescastaneus"),
set_prior("normal(0,1)", class = "sd", group = "Host", resp = "Adelgespiceae"),
set_prior("normal(0,1)", class = "sd", group = "Host", resp = "Megastigmusrafni"),
set_prior("normal(0,1)", class = "sd", group = "Host", resp = "Megastigmussuspectus"),
set_prior("normal(0,1)", class = "sd", group = "Host", resp = "Fomesfomentarius"),
set_prior("normal(0,1)", class = "sd", group = "Host", resp = "Cydiasplendana"),
set_prior("normal(0,1)", class = "sd", group = "Host", resp = "Saperdascalaris"),
set_prior("normal(0,1)", class = "sd", group = "Host", resp = "Caliroaannulipes"),
set_prior("normal(0,1)", class = "sd", group = "Host", resp = "Tortrixviridana"),
set_prior("normal(0,1)", class = "sd", group = "Host", resp = "Xiphydriacamelus"),
set_prior("normal(0,1)", class = "sd", group = "Host", resp = "Scolytusscolytus"),
set_prior("normal(0,1)", class = "sd", group = "Host", resp = "Paranthrenetabaniformis"),
set_prior("normal(0,1)", class = "sd", group = "Host", resp = "Fomitopsispinicola"),
set_prior("normal(0,1)", class = "sd", group = "Host", resp = "Xyleborusdispar"),
set_prior("normal(0,1)", class = "sd", group = "Host", resp = "Heterobasidionannosum"),
set_prior("normal(0,1)", class = "sd", group = "Host", resp = "Saturniapyri"),
set_prior("normal(0,1)", class = "sd", group = "Host", resp = "Neonectriaditissima"),
set_prior("normal(0,1)", class = "sd", group = "Host", resp = "Armillariaostoyae"),
set_prior("normal(0,1)", class = "sd", group = "Host", resp = "Hylurgopspalliatus"),
set_prior("normal(0,1)", class = "sd", group = "Host", resp = "Chionaspissalicis"),
set_prior("normal(0,1)", class = "sd", group = "Host", resp = "Heterobasidionparviporum"),
set_prior("normal(0,1)", class = "sd", group = "Host", resp = "Pityogeneschalcographus"),
set_prior("normal(0,1)", class = "sd", group = "Host", resp = "Sirococcusconigenus"),
set_prior("normal(0,1)", class = "sd", group = "Host", resp = "Gremmeniellaabietina"),
set_prior("normal(0,1)", class = "sd", group = "Host", resp = "Hylastesater"),
set_prior("normal(0,1)", class = "sd", group = "Host", resp = "Mycosphaerellapini"),
set_prior("normal(0,1)", class = "sd", group = "Host", resp = "Sphaeropsissapinea"),
set_prior("normal(0,1)", class = "sd", group = "Host", resp = "Pissodescastaneus"))
 ## Train model, using full dataset with NAs included in HSM
       fullfit <- brm(fullform, data = host_matrix, data2 = list(A=A) , prior = priors,
                      chains = 1, init = 0,
                      iter=1000, warmup=500,
                      backend = "cmdstanr",
                      threads = threading(8)
       )
       
       ## Generate predictions for the new values
       logit_df <- as.data.frame(posterior_summary(posterior_linpred(fullfit, transform = FALSE)))
    predictions_NNT<- cbind(host_matrix_NA$Host,logit_df[,seq(1,length(logit_df),4)])
    lowerquant_NNT<- cbind(host_matrix_NA$Host,logit_df[,seq(3,length(logit_df),4)])
    
    colnames(predictions_NNT)[1]<-"Host"
    colnames(predictions_NNT)[-1]<- cluster
 # Save predictions
write.csv(predictions_NNT, "~/Outputs/predictions_Adelges_piceae.csv")

