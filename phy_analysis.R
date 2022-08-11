# Thu Aug 04 10:52:59 2022 ------------------------------
# Anurans islands biogeography: how, who and why? 

# loading packages
#install.packages("adiv")
library(adiv) # package of measures biodiversity 
library(tidyverse) #data handling
library(phytools) # phylogenetic data handling 
library(picante) #checking phylogeny name
library(visdat) #check missing data

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

# Herein, we preparing data 
# community data 
dir()
anura_islands_orig <- as_tibble(read.csv2("comunidades_ilhas.csv", sep=","))
head(anura_islands_orig)
anura_islands <- anura_islands_orig %>% select(-X, -OBJECTID_1) #copy

#check data
head(anura_islands_orig) # 2505 columns with island id  
head(anura_islands) # 2503 spp without island id

# we use copy data without island id columns 
nrow(anura_islands) # 2289 islands 
ncol(anura_islands) # 2503 anuran islands spp.

# traits data
anura_traits <- read.csv2("functional_anura.csv", sep=",") 
names(anura_traits) #2249
dim(anura_traits)

# selecting traits
traits <- c("Species", "Ter","Fos","Aqu","Arb", "Dir", "Lar", "Viv", "Litter_size_max_n")
anura_traits_select <- anura_traits %>% select(traits)

# Thu Aug 11 12:05:59 2022 ------------------------------
# Join trait data with community data

#Conferindo quais espécies da lista possuem traits
# Lista para join
# Arrumando para linhas
anuran_list <- anura_islands %>%
  gather(key="Species",value="count", Abavorana_luctuosa:Zhangixalus_yinggelingensis) %>% distinct(Species, .keep_all = TRUE)

nrow(anuran_list) #2503 spp composição de espécies
nrow(anura_traits_select) # 2249 spp traits

# Left join: espécies que estão no x mantidas e agrupa com os correspiondentes de y 
anura_with_traits <- left_join(anuran_list, anura_traits_select, by="Species") %>% select(-count)
vis_miss(anura_with_traits)
# habitat use (ter, fos, aqu, arb) multi-choice traits with na's
# breeding strategy (dir, lar, viv) choice trait - 10.89% missing data
# litter_size_max continuous trait - 68% missing data.  
# View(anura_traits_select)

# Saving traits data
write.csv2(anura_with_traits, "anura_islands_traits.csv", sep=",")
?write.csv

# Anti join: espécies que estão no x, mas não estão no y (x,y)
anura_without_traits <- anti_join(anuran_list, anura_traits_select, by="Species")
nrow(anura_without_traits) # 255 spp.

#Retirando espécies da comunidade que não estão nos traits
rem.col.list <- anura_without_traits$Species
short_anura_islands <- anura_islands[,!(names(anura_islands)%in% rem.col.list)]
ncol(short_anura_islands) #2248 species with traits

###---
# phylogenetic data 
amphibia_phy <- read.tree("amph_shl_new_Consensus_7238.tre")

# we changed names in phylogeny for match to traits and composition matrix
name_phy <- as.tibble(amphibia_phy$tip.label) %>% dplyr::mutate(new_name = amphibia_phy$tip.label)
View(name_phy)

# nomenclature adjustments "Boana_semilineatus"

name_phy$new_name <- name_phy$new_name %>% str_replace_all("Hypsiboas", "Boana") %>% str_replace_all("Boana_albomarginatus", "Boana_albomarginata") %>% str_replace_all("Boana_calcaratus", "Boana_calcarata") %>% str_replace_all("Boana_fasciatus", "Boana_fasciata") %>% str_replace_all("Boana_geographicus", "Boana_geographica") %>% str_replace_all("Boana_pulchellus", "Boana_pulchella") %>% str_replace_all("Boana_multifasciatus", "Boana_multifasciata") %>% str_replace_all("Boana_ornatissimus", "Boana_ornatissima") %>% str_replace_all("Boana_punctatus", "Boana_punctata") %>% str_replace_all("Boana_rufitelus", "Boana_rufitela") %>% str_replace_all("Boana_semilineatus","Boana_semilineata") %>% str_replace_all("Boana_semilineatus", "Boana_semilineata") %>% str_replace_all("Chiromantis_doriae", "Chirixalus_doriae") %>% str_replace_all("Hylomantis_aspera", "Agalychnis_aspera") %>% str_replace_all("Pachymedusa_dacnicolor", "Agalychnis_dacnicolor") %>% str_replace_all("Hylarana_luctuosa", "Abavorana_luctuosa") %>% str_replace_all("Ingerana_baluensis", "Alcalus_baluensis") %>% str_replace_all("Ingerana_mariae", "Alcalus_mariae") %>% str_replace_all("Ingerana_rajae","Alcalus_rajae") %>% str_replace_all("Ingerana_sariba","Alcalus_sariba") %>% str_replace_all("Stumpffia_helenae","Anilany_helenae") %>% str_replace_all("Cophyla_phyllodactyla","Platypelis_phyllodactyla") %>% str_replace_all("Phyllomedusa_tomopterna", "Callimedusa_tomopterna") %>% str_replace_all("Syncope_hudsoni", "Chiasmocleis_hudsoni") %>% str_replace_all("Cyclorana_alboguttata", "Litoria_alboguttata") %>% str_replace_all("Cyclorana_australis", "Litoria_australis") %>% str_replace_all("Cyclorana_brevipes", "Litoria_brevipes")  %>% str_replace_all("Cyclorana_cryptotis", "Litoria_cryptotis") %>% str_replace_all("Cyclorana_longipes", "Litoria_longipes") %>% str_replace_all("Cyclorana_novaehollandiae", "Litoria_novaehollandiae") %>% str_replace_all("Cyclorana_nullicedens", "Litoria_nullicedens") %>% str_replace_all("Cyclorana_vagitus", "Litoria_vagitus") %>% str_replace_all("Scinax_trapicheiroi", "Ololygon_trapicheiroi") %>% str_replace_all("Scinax_agilis", "Ololygon_agilis") %>% str_replace_all("Scinax_argyreornatus", "Ololygon_argyreornata") %>% str_replace_all("Scinax_berthae", "Ololygon_berthae") %>% str_replace_all("Scinax_brieni", "Ololygon_brieni") %>% str_replace_all("Scinax_catharinae", "Ololygon_catharinae") %>% str_replace_all("Scinax_flavoguttatus", "Ololygon_flavoguttata") %>% str_replace_all("Scinax_humilis", "Ololygon_humilis") %>% str_replace_all("Scinax_jureia", "Ololygon_jureia") %>% str_replace_all("Scinax_littoralis", "Ololygon_littoralis") %>% str_replace_all("Scinax_perpusillus", "Ololygon_perpusilla") %>% str_replace_all("Scinax_rizibilis", "Ololygon_rizibilis") %>% str_replace_all("Rana_palmipes", "Lithobates_palmipes") %>% str_replace_all("Rana_berlandieri", "Lithobates_berlandieri") %>% str_replace_all("Rana_capito", "Lithobates_capito") %>% str_replace_all("Rana_catesbeiana", "Lithobates_catesbeianus") %>% str_replace_all("Rana_clamitans", "Lithobates_clamitans") %>% str_replace_all("Rana_forreri", "Lithobates_forreri") %>% str_replace_all("Rana_grylio", "Lithobates_grylio")  %>% str_replace_all("Rana_heckscheri", "Lithobates_heckscheri") %>% str_replace_all("Rana_magnaocularis", "Lithobates_magnaocularis")  %>% str_replace_all("Rana_palustris", "Lithobates_palustris")  %>% str_replace_all("Rana_pipiens", "Lithobates_pipiens")  %>% str_replace_all("Rana_septentrionalis", "Lithobates_septentrionalis") %>% str_replace_all("Rana_sphenocephala", "Lithobates_sphenocephalus") %>% str_replace_all("Rana_sylvatica", "Lithobates_sylvaticus") %>% str_replace_all("Rana_vaillanti", "Lithobates_vaillanti") %>% str_replace_all("Rana_virgatipes", "Lithobates_virgatipes")  %>% str_replace_all("Rana_warszewitschii", "Lithobates_warszewitschii") %>% str_replace_all("Rana_yavapaiensis", "Lithobates_yavapaiensis") %>% str_replace_all("Lysapsus_boliviana", "Lysapsus_bolivianus") 

name_phy$new_name <- name_phy$new_name %>% str_replace_all("Fejervarya_andamanensis", "Minervarya_andamanensis") %>% str_replace_all("Fejervarya_greenii", "Minervarya_greenii") %>% str_replace_all("Fejervarya_nepalensis", "Minervarya_nepalensis") %>% str_replace_all("Fejervarya_nicobariensis", "Minervarya_nicobariensis") %>% str_replace_all("Fejervarya_pierrei", "Minervarya_pierrei") %>% str_replace_all("Fejervarya_rufescens", "Minervarya_rufescens") %>% str_replace_all("Fejervarya_teraiensis", "Minervarya_teraiensis") %>% str_replace_all("Leptolalax_arayai", "Leptobrachella_arayai") %>% str_replace_all("Leptolalax_dringi", "Leptobrachella_dringi") %>% str_replace_all("Leptolalax_fuliginosus", "Leptobrachella_fuliginosa" ) %>% str_replace_all("Leptolalax_gracilis","Leptobrachella_gracilis") %>% str_replace_all("Leptolalax_hamidi","Leptobrachella_hamidi") %>% str_replace_all("Leptolalax_kajangensis","Leptobrachella_kajangensis") %>% str_replace_all("Leptolalax_liui","Leptobrachella_liui") %>% str_replace_all("Leptolalax_maurus", "Leptobrachella_maura") %>% str_replace_all("Leptolalax_pelodytoides","Leptobrachella_pelodytoides") %>% str_replace_all("Leptolalax_picta","Leptobrachella_pictus") %>% str_replace_all("Fejervarya_kirtisinghei", "Minervarya_kirtisinghei") %>% str_replace_all("Liophryne_allisoni","Sphenophryne_allisoni") %>% str_replace_all("Liophryne_dentata", "Sphenophryne_dentata") %>% str_replace_all("Liophryne_magnitympanum","Sphenophryne_magnitympanum") %>% str_replace_all("Liophryne_rhododactyla","Sphenophryne_rhododactyla") %>% str_replace_all("Liophryne_rubra","Sphenophryne_rubra") %>% str_replace_all("Liophryne_rubra","Sphenophryne_rubra") %>% str_replace_all("Liophryne_schlaginhaufeni", "Sphenophryne_schlaginhaufeni") %>% str_replace_all("Liophryne_similis","Sphenophryne_similis") %>% str_replace_all("Oxydactyla_brevicrus","Sphenophryne_brevicrus") %>% str_replace_all("Oxydactyla_coggeri","Sphenophryne_coggeri") %>% str_replace_all("Oxydactyla_crassa","Sphenophryne_crassa") %>% str_replace_all("Oxydactyla_stenodactyla","Sphenophryne_stenodactyla") %>% str_replace_all("Rhacophorus_achantharrhena", "Zhangixalus_achantharrhena") %>% str_replace_all("Rhacophorus_arboreus", "Zhangixalus_arboreus") %>% str_replace_all("Rhacophorus_arvalis", "Zhangixalus_arvalis") %>% str_replace_all("Rhacophorus_aurantiventris", "Zhangixalus_aurantiventris") %>% str_replace_all("Rhacophorus_dennysi", "Zhangixalus_dennysi") %>% str_replace_all("Rhacophorus_dulitensis", "Zhangixalus_dulitensis") %>% str_replace_all("Rhacophorus_moltrechti", "Zhangixalus_moltrechti") %>% str_replace_all("Rhacophorus_owstoni", "Zhangixalus_owstoni") %>% str_replace_all("Rhacophorus_prasinatus", "Zhangixalus_prasinatus") %>% str_replace_all("Rhacophorus_prominanus", "Zhangixalus_prominanus") %>% str_replace_all("Rhacophorus_schlegelii", "Zhangixalus_schlegelii") %>% str_replace_all("Rhacophorus_taipeianus", "Zhangixalus_taipeianus") %>% str_replace_all("Rhacophorus_viridis", "Zhangixalus_viridis") %>% str_replace_all("Platymantis", "Cornufer")

amphibia_phy$tip.label <- name_phy$new_name

# Checking species match
options(max.print=250)
match.phylo.comm(amphibia_phy,short_anura_islands)

# genus
names(anura_islands) <- gsub("Amnirana", "Hylarana", fixed=T, names(anura_islands)) 
names(anura_islands) <- gsub("Cophyla", "Platypelis", fixed=T, names(anura_islands)) 
names(anura_islands) <- gsub("Dryophytes", "Hyla", fixed=T, names(anura_islands)) 
names(anura_islands) <- gsub("Chalcorana", "Hylarana", fixed=T, names(anura_islands)) 

#Prune with phylogeny Pyron and Wiens consensus
phy_anura_islands <- prune.sample(short_anura_islands, amphibia_phy) 
# 2159 spp
phy_anura_islands

#quantas espécies perdemos?
2248-2041 # 207 # ncol(short_anura_islands) #2248

# Phylogeny anuran islands
plot(phy_anura, type="fan",show.tip.label = FALSE)
#save phy islands
write.tree(phy_anura, file = "anura_islands.tre")

###
[1] "Dropping taxa from the community because they are not present in the phylogeny:"
[1] "Amnirana_albolabris"          "Amnirana_amnicola"           
[3] "Amnirana_galamensis"          "Amnirana_lepus"              
[5] "Amnirana_occidentalis"        "Anodonthyla_boulengeri"      
[7] "Aphantophryne_nana"           "Asterophrys_eurydactyla"     
[9] "Asterophrys_marani"           "Asterophrys_pullifer"        
[11] "Asterophrys_slateri"          "Bijurana_nicobariensis"      

            "Chalcorana_chalconota"       
[15] "Chalcorana_crassiovis"        "Chalcorana_kampeni"          
[17] "Chalcorana_macrops"           "Chalcorana_megalonesa"       
[19] "Chalcorana_mocquardii"        "Chalcorana_parvaccola"       
[21] "Chalcorana_raniceps"          "Chalcorana_rufipes"          

            [23] "Chirixalus_nongkhorensis"     

"Choerophryne_brunhildae"     
[25] "Choerophryne_darlingtoni"     "Choerophryne_exclamitans"    
[27] "Choerophryne_fafniri"         "Choerophryne_gunnari"        
[29] "Choerophryne_laurini"         "Choerophryne_rhenaurum"      
[31] "Choerophryne_sanguinopicta"   "Choerophryne_siegfriedi"     
[33] "Choerophryne_swanhildae"      "Choerophryne_tubercula"      
[35] "Choerophryne_valkuriarum"     "Choerophryne_variegata"      

[37] "Cophyla_barbouri"             "Cophyla_grandis"             
[39] "Cophyla_mavomavo"             "Cophyla_milloti"             
[41] "Cophyla_phyllodactyla"        "Cophyla_pollicaris"          
[43] "Cophyla_tetra"                "Cophyla_tsaratananaensis"    
[45] "Cophyla_tuberifera"           

"Copiula_alpestris"           
[47] "Copiula_guttata"              

"Cornufer_elegans"            
[49] "Cornufer_gigas"               "Cornufer_hedigeri"           
[51] "Cornufer_heffernani"          "Cornufer_malukuna"           
[53] "Cornufer_mediodiscus"         "Cornufer_minutus"            
[55] "Cornufer_opisthodon"          "Cornufer_trossulus"          
[57] "Cornufer_vertebralis"         "Cornufer_vogti"              
[59] "Cornufer_wolfi"               

"Dryophytes_andersonii"       
[61] "Dryophytes_avivoca"           "Dryophytes_chrysoscelis"     
[63] "Dryophytes_cinereus"          "Dryophytes_femoralis"        
[65] "Dryophytes_gratiosus"         "Dryophytes_immaculatus"      
[67] "Dryophytes_japonicus"         "Dryophytes_squirellus"       
[69] "Dryophytes_versicolor"        

"Fejervarya_pulla"            
[71] "Glandirana_emeljanovi"        "Glandirana_rugosa"           

[73] "Glyphoglossus_brooksii"       "Glyphoglossus_flavus"        
[75] "Glyphoglossus_guttulatus"     "Glyphoglossus_smithi"        
[77] "Glyphoglossus_volzi"          "Hydrophylax_gracilis"        
[79] "Hydrophylax_leptoglossa"      "Hydrophylax_malabaricus"     
[81] "Hylophorbus_infulatus"        "Hyperolius_koehleri"         
[83] "Indosylvirana_aurantiaca"     "Indosylvirana_temporalis"    
[85] "Leptobrachella_picta"         "Limnonectes_rhacodus"        
[87] "Litoria_rueppelli"            "Mantophryne_menziesi"        

[89] "Megophrys_aceras"             "Megophrys_baluensis"         
[91] "Megophrys_boettgeri"          "Megophrys_brachykolos"       
[93] "Megophrys_carinense"          "Megophrys_dringi"            
[95] "Megophrys_edwardinae"         "Megophrys_major"             
[97] "Megophrys_parallela"          "Megophrys_parva"             

[99] "Micryletta_erythropoda"       "Minervarya_syhadrensis"      
[101] "Nanohyla_perparva"            "Nanohyla_petrigena"          
[103] "Nidirana_adenopleura"         "Nidirana_hainanensis"        
[105] "Nidirana_okinavana"           "Nyctimystes_infrafrenatus"   
[107] "Nyctimystes_sauroni"          "Nyctimystes_tyleri"          

[109] "Papurana_arfaki"              "Papurana_attigua"            
[111] "Papurana_aurata"              "Papurana_daemeli"            
[113] "Papurana_elberti"             "Papurana_florensis"          
[115] "Papurana_garritor"            "Papurana_grisea"             
[117] "Papurana_jimiensis"           "Papurana_kreffti"            
[119] "Papurana_milleti"             "Papurana_milneana"           
[121] "Papurana_moluccana"           "Papurana_novaeguineae"       
[123] "Papurana_papua"               "Papurana_supragrisea"        
[125] "Papurana_volkerjane"          "Papurana_waliesa"            
[127] "Phrynoidis_asper"             "Pithecopus_hypochondrialis"  

[129] "Platymantis_banahao"          "Platymantis_bayani"          
[131] "Platymantis_biak"             "Platymantis_cagayanensis"    
[133] "Platymantis_cornutus"         "Platymantis_corrugatus"      
[135] "Platymantis_diesmosi"         "Platymantis_dorsalis"        
[137] "Platymantis_guentheri"        "Platymantis_hazelae"         
[139] "Platymantis_indeprensus"      "Platymantis_isarog"          
[141] "Platymantis_lawtoni"          "Platymantis_levigatus"       
[143] "Platymantis_luzonensis"       "Platymantis_mimulus"         
[145] "Platymantis_montanus"         "Platymantis_naomii"          
[147] "Platymantis_negrosensis"      "Platymantis_paengi"          
[149] "Platymantis_panayensis"       "Platymantis_polillensis"     
[151] "Platymantis_pseudodorsalis"   "Platymantis_pygmaeus"        
[153] "Platymantis_rabori"           "Platymantis_sierramadrensis" 
[155] "Platymantis_spelaeus"         "Platymantis_subterrestris"   
[157] "Platymantis_taylori"          

"Plethodontohyla_alluaudi"    

[159] "Pulchrana_baramica"           "Pulchrana_debussyi"          
[161] "Pulchrana_glandulosa"         "Pulchrana_grandocula"        
[163] "Pulchrana_laterimaculata"     "Pulchrana_mangyanum"         
[165] "Pulchrana_melanomenta"        "Pulchrana_moellendorffi"     
[167] "Pulchrana_picturata"          "Pulchrana_siberu"            
[169] "Pulchrana_signata"            "Pulchrana_similis"           

[171] "Rentapia_everetti"            "Rentapia_hosii"              
[173] "Rhacophorus_rhyssocephalus"   "Rohanixalus_hansenae"        
[175] "Rohanixalus_vittatus"         "Sanguirana_mearnsi"          

[177] "Sclerophrys_buchneri"         "Sclerophrys_camerunensis"    
[179] "Sclerophrys_funerea"          "Sclerophrys_gracilipes"      
[181] "Sclerophrys_gutturalis"       "Sclerophrys_kassasii"        
[183] "Sclerophrys_latifrons"        "Sclerophrys_maculata"        
[185] "Sclerophrys_pentoni"          "Sclerophrys_regularis"       
[187] "Sclerophrys_superciliaris"    "Sclerophrys_tihamica"        
[189] "Sclerophrys_togoensis"        "Sclerophrys_tuberosa"        
[191] "Sclerophrys_xeros"            

[203] "Sphenophryne_thomsoni" 

"Sylvirana_cubitalis"         
[205] "Sylvirana_guentheri"          "Sylvirana_maosonensis"       
[207] "Sylvirana_mortenseni"         "Sylvirana_nigrovittata"      
[209] "Sylvirana_spinulosa"          

"Trachycephalus_typhonius" 

[211] "Uperodon_nagaoi"              "Uperodon_obscurus"           
[213] "Uperodon_palmatus"            "Uperodon_taprobanicus"       
[215] "Uperodon_variegatus"          

"Wijayarana_masonii"          
[217] "Wijayarana_melasma"           "Wijayarana_modiglianii"      
[219] "Wijayarana_sumatrana"         
