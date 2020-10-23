
dlbc_maf <- GDCquery_Maf("DLBC", pipelines = "mutect")

data.frame("Missense Variant" = 
             str_detect(dlbc_maf$all_effects, "missense_variant")) 


