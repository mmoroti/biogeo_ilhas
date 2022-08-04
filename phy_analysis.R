# Thu Aug 04 10:52:59 2022 ------------------------------
# Anurans islands biogeography: how, who and why? 

# loading packages
#install.packages("adiv")
library(adiv) # package of measures biodiversity 
library(tidyverse) #data handling
library(phytools) # phylogenetic data handling 
library(picante) #checking phylogeny name

# exploring 'adiv' package
data(batcomm)
batcomm$ab
batcomm$ab2
batcomm$tre

data(batcomm)
abgdivparam(batcomm$ab)
plot(abgdivparam(batcomm$ab))

abgdivparam(batcomm$ab, q=0:4)
plot(abgdivparam(batcomm$ab, q=0:4))
View(batcomm$ab)

# community data 
# how to create
dir()
anura_islands <- read.csv2("comunidades_ilhas.csv", sep=",")
nrow(anura_islands) # 2289 islands
ncol(anura_islands) # 2203 anuran islands spp.

# phylogenetic data 
amphibia_phy <- read.tree("amph_shl_new_Consensus_7238.tre")
name_phy <- as.data.frame(amphibia_phy$tip.label)
View(name_phy)

# Checking species match
match.phylo.comm(amphibia_phy,anura_islands)
options(max.print=100)

# nomenclature adjustments 
# genus
names(anura_islands) <- gsub("Hylomantis","Agalychnis", fixed=T, names(anura_islands)) 
names(anura_islands) <- gsub("Boana", "Hypsiboas", fixed=T, names(anura_islands)) 

# epitete
names(anura_islands) <- gsub("Hypsiboas_albomarginata", "Hypsiboas_albomarginatus", fixed=T, names(anura_islands)) 
names(anura_islands) <- gsub("Hypsiboas_calcarata", "Hypsiboas_calcaratus", fixed=T, names(anura_islands)) 
names(anura_islands) <- gsub("Hypsiboas_fasciata", "Hypsiboas_fasciatus", fixed=T, names(anura_islands)) 
names(anura_islands) <- gsub("Hypsiboas_geographica", "Hypsiboas_greographicus", fixed=T, names(anura_islands)) 
names(anura_islands) <- gsub("Hypsiboas_multifasciata", "Hypsiboas_multifasciatus", fixed=T, names(anura_islands)) 
names(anura_islands) <- gsub("Hypsiboas_ornatissima", "Hypsiboas_ornatissimus", fixed=T, names(anura_islands)) 
names(anura_islands) <- gsub("Hypsiboas_punctata", "Hypsiboas_punctatus", fixed=T, names(anura_islands)) 
names(anura_islands) <- gsub("Hypsiboas_rufitela", "Hypsiboas_rufitelus", fixed=T, names(anura_islands)) 
names(anura_islands) <- gsub("Hypsiboas_semilineata", "Hypsiboas_semilineatus", fixed=T, names(anura_islands)) 
names(anura_islands) <- gsub("Chirixalus_doriae", "Chiromantis_doriae", fixed=T, names(anura_islands)) 

names(anura_islands) <- gsub("Agalychnis_aspera", "Hylomantis_aspera", fixed=T, names(anura_islands))



# Prune with phylogeny Pyron and Wiens consensus
phy_anura <- prune.sample(anura_islands, amphibia_phy) # 2029 spp
phy_anura


# spp without phylogeny 
# "Choerophryne_alainduboisi", "Choerophryne_alpestris","Choerophryne_bickfordi","Choerophryne_bisyllaba","Choerophryne_brevicrus","Choerophryne_brunhildae","Choerophryne_crucifer","Choerophryne_darlingtoni","Choerophryne_epirrhina","Choerophryne_exclamitans","Choerophryne_fafniri","Choerophryne_grylloides","Choerophryne_gudrunae","Choerophryne_gunnari","Choerophryne_laurini","Choerophryne_multisyllaba","Choerophryne_murrita","Choerophryne_pandanicola","Choerophryne_pipiens","Choerophryne_rhenaurum","Choerophryne_sanguinopicta","Choerophryne_siegfriedi","Choerophryne_swanhildae","Choerophryne_tubercula","Choerophryne_valkuriarum","Choerophryne_variegata"


[1] "Dropping taxa from the community because they are not present in the phylogeny:"
[1] "X"                            "OBJECTID_1"                  
[3] "Abavorana_luctuosa"           "Agalychnis_aspera"           
[5] "Agalychnis_dacnicolor"        "Aglyptodactylus_australis"   
[7] "Aglyptodactylus_chorus"       "Aglyptodactylus_inguinalis"  
[9] "Alcalus_baluensis"            "Alcalus_mariae"              
[11] "Alcalus_rajae"                "Alcalus_sariba"              
[13] "Amnirana_albolabris"          "Amnirana_amnicola"           
[15] "Amnirana_galamensis"          "Amnirana_lepus"              
[17] "Amnirana_occidentalis"        "Amolops_albispinus"          
[19] "Anilany_helenae"              "Anodonthyla_boulengeri"      
[21] "Ansonia_phuketensis"          "Ansonia_teneritas"           
[23] "Ansonia_vidua"                "Aphantophryne_nana"          
[25] "Asterophrys_eurydactyla"      "Asterophrys_foja"            
[27] "Asterophrys_marani"           "Asterophrys_pullifer"        
[29] "Asterophrys_slateri"          "Austrochaperina_alexanderi"  
[31] "Austrochaperina_laurae"       "Austrochaperina_punctata"    
[33] "Austrochaperina_rudolfarndti" "Bijurana_nicobariensis"      
[35] "Blommersia_nataliae"          "Blommersia_transmarina"      
[37] "Hypsiboas_greographicus"      "Hypsiboas_pulchella"         
[39] "Boophis_ankarafensis"         "Boophis_boppa"               
[41] "Boophis_nauticus"             "Buergeria_choui"             
[43] "Buergeria_otai"               "Callimedusa_tomopterna"      
[45] "Callulops_argus"              "Callulops_biakensis"         
[47] "Callulops_bicolor"            "Callulops_neuhaussi"         
[49] "Callulops_stellatus"          "Callulops_taxispilotus"      
[51] "Callulops_wondiwoiensis"      "Callulops_yapenensis"        
[53] "Chalcorana_chalconota"        "Chalcorana_crassiovis"       
[55] "Chalcorana_kampeni"           "Chalcorana_macrops"          
[57] "Chalcorana_megalonesa"        "Chalcorana_mocquardii"       
[59] "Chalcorana_parvaccola"        "Chalcorana_raniceps"         
[61] "Chalcorana_rufipes"           "Chiasmocleis_hudsoni"        
[63] "Chirixalus_nongkhorensis"     "Chirixalus_trilaksonoi"      

[91] "Cophixalus_cateae"            "Cophixalus_hannahae"         
[93] "Cophixalus_rajampatensis"     "Cophixalus_salawatiensis"    
[95] "Cophixalus_viridis"           "Cophyla_alticola"            
[97] "Cophyla_barbouri"             "Cophyla_grandis"             
[99] "Cophyla_karenae"              "Cophyla_maharipeo" 