# wxGenerateAnatomicalFeatures
# This function creates a list with indexed morphological features
# José Luis Alatorre Warren, University of Zurich, 2017

wxGenerateAnatomicalFeatures <- function(mainDataset) {
  
  # Create a list with all anatomical features
  anatomicalFeatures <- list()
  
  # Create one array per feature (with all points of all specimens)
  anatomicalFeatures[["hpb_cerebellum_delimitation_left"]]$dense             <- mainDataset$rotated[   1:  98,,]
  anatomicalFeatures[["hpb_cerebellum_delimitation_right"]]$dense            <- mainDataset$rotated[  99: 196,,]
  anatomicalFeatures[["hpb_fissura_longitudinalis_cerebri_frontal"]]$dense   <- mainDataset$rotated[ 197: 335,,]
  anatomicalFeatures[["hpb_fissura_longitudinalis_cerebri_occipital"]]$dense <- mainDataset$rotated[ 336: 383,,]
  anatomicalFeatures[["hpb_fissura_longitudinalis_cerebri_parietal"]]$dense  <- mainDataset$rotated[ 384: 423,,]
  anatomicalFeatures[["hpb_sulcus_centralis_left"]]$dense                    <- mainDataset$rotated[ 424: 508,,]
  anatomicalFeatures[["hpb_sulcus_centralis_right"]]$dense                   <- mainDataset$rotated[ 509: 593,,]
  anatomicalFeatures[["hpb_sulcus_frontalis_inferior_left"]]$dense           <- mainDataset$rotated[ 594: 628,,]
  anatomicalFeatures[["hpb_sulcus_frontalis_inferior_right"]]$dense          <- mainDataset$rotated[ 629: 663,,]
  anatomicalFeatures[["hpb_sulcus_frontalis_superior_left"]]$dense           <- mainDataset$rotated[ 664: 748,,]
  anatomicalFeatures[["hpb_sulcus_frontalis_superior_right"]]$dense          <- mainDataset$rotated[ 749: 833,,]
  anatomicalFeatures[["hpb_sulcus_intraparietalis_left"]]$dense              <- mainDataset$rotated[ 834: 872,,]
  anatomicalFeatures[["hpb_sulcus_intraparietalis_right"]]$dense             <- mainDataset$rotated[ 873: 911,,]
  anatomicalFeatures[["hpb_sulcus_lateralis_left"]]$dense                    <- mainDataset$rotated[ 912:1014,,]
  anatomicalFeatures[["hpb_sulcus_lateralis_right"]]$dense                   <- mainDataset$rotated[1015:1117,,]
  anatomicalFeatures[["hpb_sulcus_postcentralis_left"]]$dense                <- mainDataset$rotated[1118:1178,,]
  anatomicalFeatures[["hpb_sulcus_postcentralis_right"]]$dense               <- mainDataset$rotated[1179:1239,,]
  anatomicalFeatures[["hpb_sulcus_precentralis_left"]]$dense                 <- mainDataset$rotated[1240:1316,,]
  anatomicalFeatures[["hpb_sulcus_precentralis_right"]]$dense                <- mainDataset$rotated[1317:1393,,]
  anatomicalFeatures[["hpb_sulcus_temporalis_superior_left"]]$dense          <- mainDataset$rotated[1394:1483,,]
  anatomicalFeatures[["hpb_sulcus_temporalis_superior_right"]]$dense         <- mainDataset$rotated[1484:1573,,]
  anatomicalFeatures[["hpi_alae_parvae_left"]]$dense                         <- mainDataset$rotated[1574:1605,,]
  anatomicalFeatures[["hpi_alae_parvae_right"]]$dense                        <- mainDataset$rotated[1606:1637,,]
  anatomicalFeatures[["hpi_basion_dorsum_sella"]]$dense                      <- mainDataset$rotated[1638:1676,,]
  anatomicalFeatures[["hpi_crista_occipitalis_interna"]]$dense               <- mainDataset$rotated[1677:1711,,]
  anatomicalFeatures[["hpi_foramen_magnum"]]$dense                           <- mainDataset$rotated[1712:1807,,]
  anatomicalFeatures[["hpi_margo_superior_partis_petrosae_left"]]$dense      <- mainDataset$rotated[1808:1861,,]
  anatomicalFeatures[["hpi_margo_superior_partis_petrosae_right"]]$dense     <- mainDataset$rotated[1862:1915,,]
  anatomicalFeatures[["hpi_midline_parietal_bone"]]$dense                    <- mainDataset$rotated[1916:2014,,]
  anatomicalFeatures[["hpi_midline_frontal_bone"]]$dense                     <- mainDataset$rotated[2015:2116,,]
  anatomicalFeatures[["hpi_midline_occipital_bone"]]$dense                   <- mainDataset$rotated[2117:2162,,]
  anatomicalFeatures[["hpi_sella_turcica_foramen_caecum"]]$dense             <- mainDataset$rotated[2163:2210,,]
  anatomicalFeatures[["hpi_sinus_transversus_sigmoideus_left"]]$dense        <- mainDataset$rotated[2211:2310,,]
  anatomicalFeatures[["hpi_sinus_transversus_sigmoideus_right"]]$dense       <- mainDataset$rotated[2311:2410,,]
  anatomicalFeatures[["hpo_asterion_pterion_left"]]$dense                    <- mainDataset$rotated[2411:2503,,]
  anatomicalFeatures[["hpo_asterion_pterion_right"]]$dense                   <- mainDataset$rotated[2504:2596,,]
  anatomicalFeatures[["hpo_bregma_nasion"]]$dense                            <- mainDataset$rotated[2597:2704,,]
  anatomicalFeatures[["hpo_lambda_opisthion"]]$dense                         <- mainDataset$rotated[2705:2797,,]
  anatomicalFeatures[["hpo_sutura_coronalis_left"]]$dense                    <- mainDataset$rotated[2798:2897,,]
  anatomicalFeatures[["hpo_sutura_coronalis_right"]]$dense                   <- mainDataset$rotated[2898:2997,,]
  anatomicalFeatures[["hpo_sutura_lambdoidea_left"]]$dense                   <- mainDataset$rotated[2998:3075,,]
  anatomicalFeatures[["hpo_sutura_lambdoidea_right"]]$dense                  <- mainDataset$rotated[3076:3153,,]
  anatomicalFeatures[["hpo_sutura_sagittalis"]]$dense                        <- mainDataset$rotated[3154:3258,,]
  
  # Resample features such that each feature has its two extremes landmarks plus every 10th landmark in between
  totalNumberOfFeatures <- length(anatomicalFeatures)
  for (ii in 1:totalNumberOfFeatures){
    anatomicalFeatures[[ii]][["sparse"]]   <- wxSimplify(anatomicalFeatures[[ii]][["dense"]],"sparse")
    anatomicalFeatures[[ii]][["extremes"]] <- wxSimplify(anatomicalFeatures[[ii]][["dense"]],"extremes")    
  }
  
  return(anatomicalFeatures)
  
}