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

# nomenclature adjustments 
name_phy$new_name <- name_phy$new_name %>% str_replace_all("Hypsiboas", "Boana") %>% str_replace_all("Boana_albomarginatus", "Boana_albomarginata") %>% str_replace_all("Boana_calcaratus", "Boana_calcarata") %>% str_replace_all("Boana_fasciatus", "Boana_fasciata") %>% str_replace_all("Boana_geographicus", "Boana_geographica") %>% str_replace_all("Boana_pulchellus", "Boana_pulchella") %>% str_replace_all("Boana_multifasciatus", "Boana_multifasciata") %>% str_replace_all("Boana_ornatissimus", "Boana_ornatissima") %>% str_replace_all("Boana_punctatus", "Boana_punctata") %>% str_replace_all("Boana_rufitelus", "Boana_rufitela") %>% str_replace_all("Boana_semilineatus", "Boana_selimineata") %>% str_replace_all("Boana_semilineatus", "Boana_semilineata") %>% str_replace_all("Chiromantis_doriae", "Chirixalus_doriae") %>% str_replace_all("Hylomantis_aspera", "Agalychnis_aspera") %>% str_replace_all("Pachymedusa_dacnicolor", "Agalychnis_dacnicolor") %>% str_replace_all("Hylarana_luctuosa", "Abavorana_luctuosa") %>% str_replace_all("Ingerana_baluensis", "Alcalus_baluensis") %>% str_replace_all("Ingerana_mariae", "Alcalus_mariae") %>% str_replace_all("Ingerana_rajae","Alcalus_rajae") %>% str_replace_all("Ingerana_sariba","Alcalus_sariba") %>% str_replace_all("Stumpffia_helenae","Anilany_helenae") %>% str_replace_all("Cophyla_phyllodactyla","Platypelis_phyllodactyla") %>% str_replace_all("Phyllomedusa_tomopterna", "Callimedusa_tomopterna") %>% str_replace_all("Syncope_hudsoni", "Chiasmocleis_hudsoni") %>% str_replace_all("Cyclorana_alboguttata", "Litoria_alboguttata") %>% str_replace_all("Cyclorana_australis", "Litoria_australis") %>% str_replace_all("Cyclorana_brevipes", "Litoria_brevipes")  %>% str_replace_all("Cyclorana_cryptotis", "Litoria_cryptotis") %>% str_replace_all("Cyclorana_longipes", "Litoria_longipes") %>% str_replace_all("Cyclorana_novaehollandiae", "Litoria_novaehollandiae") %>% str_replace_all("Cyclorana_nullicedens", "Litoria_nullicedens") %>% str_replace_all("Cyclorana_vagitus", "Litoria_vagitus") %>% str_replace_all("Scinax_trapicheiroi", "Ololygon_trapicheiroi") %>% str_replace_all("Scinax_agilis", "Ololygon_agilis") %>% str_replace_all("Scinax_argyreornatus", "Ololygon_argyreornata") %>% str_replace_all("Scinax_berthae", "Ololygon_berthae") %>% str_replace_all("Scinax_brieni", "Ololygon_brieni") %>% str_replace_all("Scinax_catharinae", "Ololygon_catharinae") %>% str_replace_all("Scinax_flavoguttatus", "Ololygon_flavoguttata") %>% str_replace_all("Scinax_humilis", "Ololygon_humilis") %>% str_replace_all("Scinax_jureia", "Ololygon_jureia") %>% str_replace_all("Scinax_littoralis", "Ololygon_littoralis") %>% str_replace_all("Scinax_perpusillus", "Ololygon_perpusilla") %>% str_replace_all("Scinax_rizibilis", "Ololygon_rizibilis") %>% str_replace_all("Rana_palmipes", "Lithobates_palmipes") %>% str_replace_all("Rana_berlandieri", "Lithobates_berlandieri") %>% str_replace_all("Rana_capito", "Lithobates_capito") %>% str_replace_all("Rana_catesbeiana", "Lithobates_catesbeianus") %>% str_replace_all("Rana_clamitans", "Lithobates_clamitans") %>% str_replace_all("Rana_forreri", "Lithobates_forreri")%>% str_replace_all("Rana_grylio", "Lithobates_grylio")  %>% str_replace_all("Rana_heckscheri", "Lithobates_heckscheri") %>% str_replace_all("Rana_magnaocularis", "Lithobates_magnaocularis")  %>% str_replace_all("Rana_palustris", "Lithobates_palustris")  %>% str_replace_all("Rana_pipiens", "Lithobates_pipiens")  %>% str_replace_all("Rana_septentrionalis", "Lithobates_septentrionalis") %>% str_replace_all("Rana_sphenocephala", "Lithobates_sphenocephalus") %>% str_replace_all("Rana_sylvatica", "Lithobates_sylvaticus") %>% str_replace_all("Rana_vaillanti", "Lithobates_vaillanti") %>% str_replace_all("Rana_virgatipes", "Lithobates_virgatipes")  %>% str_replace_all("Rana_warszewitschii", "Lithobates_warszewitschii") %>% str_replace_all("Rana_yavapaiensis", "Lithobates_yavapaiensis") %>% str_replace_all("Lysapsus_boliviana", "Lysapsus_bolivianus") %>% str_replace_all("Fejervarya_andamanensis", "Minervarya_andamanensis")  %>% str_replace_all("Fejervarya_greenii", "Minervarya_greenii") %>% str_replace_all("Fejervarya_kirtisinghei", "Minervarya_kirtisinghei") %>% str_replace_all("Fejervarya_nepalensis", "Minervarya_nepalensis") %>% str_replace_all("Fejervarya_nicobariensis", "Minervarya_nicobariensis") %>% str_replace_all("Fejervarya_pierrei", "Minervarya_pierrei") %>% str_replace_all("Fejervarya_rufescens", "Minervarya_rufescens") %>% str_replace_all("Fejervarya_teraiensis", "Minervarya_teraiensis") %>% str_replace_all("Leptolalax_arayai", "Leptobrachella_arayai") %>% str_replace_all("Leptolalax_dringi", "Leptobrachella_dringi") %>% str_replace_all("Leptolalax_fuliginosus", "Leptobrachella_fuliginosa" ) %>% str_replace_all("Leptolalax_gracilis","Leptobrachella_gracilis") %>% str_replace_all("Leptolalax_hamidi","Leptobrachella_hamidi") %>% str_replace_all("Leptolalax_kajangensis","Leptobrachella_kajangensis") %>% str_replace_all("Leptolalax_liui","Leptobrachella_liui") %>% str_replace_all("Leptolalax_maurus", "Leptobrachella_maura") %>% str_replace_all("Leptolalax_pelodytoides","Leptobrachella_pelodytoides") %>% str_replace_all("Leptolalax_picta","Leptobrachella_pictus") %>% str_replace_all( "Platymantis", "Cornufer")
  
amphibia_phy$tip.label <- name_phy$new_name

# Checking species match
options(max.print=250)
match.phylo.comm(amphibia_phy,short_anura_islands)

# genus
names(anura_islands) <- gsub("Amnirana", "Hylarana", fixed=T, names(anura_islands)) 
names(anura_islands) <- gsub("Cophyla", "Platypelis", fixed=T, names(anura_islands)) 
names(anura_islands) <- gsub("Dryophytes", "Hyla", fixed=T, names(anura_islands)) 
names(anura_islands) <- gsub("Zhangixalus", "Rhacophorus", fixed=T, names(anura_islands)) 
names(anura_islands) <- gsub("Minervarya", "Fejervarya", fixed=T, names(anura_islands)) 
names(anura_islands) <- gsub("Chalcorana", "Hylarana", fixed=T, names(anura_islands)) 
names(anura_islands) <- gsub("Leptobrachella", "Leptolalax", fixed=T, names(anura_islands))

#Prune with phylogeny Pyron and Wiens consensus
phy_anura_islands <- prune.sample(short_anura_islands, amphibia_phy) 
# 2159 spp
phy_anura_islands

#quantas espécies perdemos?
2248-2007 # 241 # ncol(short_anura_islands) #2248

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
[13] "Boana_semilineata"            "Chalcorana_chalconota"       
[15] "Chalcorana_crassiovis"        "Chalcorana_kampeni"          
[17] "Chalcorana_macrops"           "Chalcorana_megalonesa"       
[19] "Chalcorana_mocquardii"        "Chalcorana_parvaccola"       
[21] "Chalcorana_raniceps"          "Chalcorana_rufipes"          
[23] "Chirixalus_nongkhorensis"     "Choerophryne_brunhildae"     
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
[45] "Cophyla_tuberifera"           "Copiula_alpestris"           
[47] "Copiula_guttata"              "Cornufer_elegans"            
[49] "Cornufer_gigas"               "Cornufer_hedigeri"           
[51] "Cornufer_heffernani"          "Cornufer_malukuna"           
[53] "Cornufer_mediodiscus"         "Cornufer_minutus"            
[55] "Cornufer_opisthodon"          "Cornufer_trossulus"          
[57] "Cornufer_vertebralis"         "Cornufer_vogti"              
[59] "Cornufer_wolfi"               "Dryophytes_andersonii"       
[61] "Dryophytes_avivoca"           "Dryophytes_chrysoscelis"     
[63] "Dryophytes_cinereus"          "Dryophytes_femoralis"        
[65] "Dryophytes_gratiosus"         "Dryophytes_immaculatus"      
[67] "Dryophytes_japonicus"         "Dryophytes_squirellus"       
[69] "Dryophytes_versicolor"        "Fejervarya_pulla"            
[71] "Glandirana_emeljanovi"        "Glandirana_rugosa"           
[73] "Glyphoglossus_brooksii"       "Glyphoglossus_flavus"        
[75] "Glyphoglossus_guttulatus"     "Glyphoglossus_smithi"        
[77] "Glyphoglossus_volzi"          "Hydrophylax_gracilis"        
[79] "Hydrophylax_leptoglossa"      "Hydrophylax_malabaricus"     
[81] "Hylophorbus_infulatus"        "Hyperolius_koehleri"         
[83] "Indosylvirana_aurantiaca"     "Indosylvirana_temporalis"    
        
[95] "Limnonectes_rhacodus"         "Litoria_rueppelli"           
[97] "Mantophryne_menziesi"         "Megophrys_aceras"            
[99] "Megophrys_baluensis"          "Megophrys_boettgeri"         
[101] "Megophrys_brachykolos"        "Megophrys_carinense"         
[103] "Megophrys_dringi"             "Megophrys_edwardinae"        
[105] "Megophrys_major"              "Megophrys_parallela"         
[107] "Megophrys_parva"              "Micryletta_erythropoda"      
[109] "Minervarya_syhadrensis"       "Nanohyla_perparva"           
[111] "Nanohyla_petrigena"           "Nidirana_adenopleura"        
[113] "Nidirana_hainanensis"         "Nidirana_okinavana"          
[115] "Nyctimystes_infrafrenatus"    "Nyctimystes_sauroni"         
[117] "Nyctimystes_tyleri"           "Papurana_arfaki"             
[119] "Papurana_attigua"             "Papurana_aurata"             
[121] "Papurana_daemeli"             "Papurana_elberti"            
[123] "Papurana_florensis"           "Papurana_garritor"           
[125] "Papurana_grisea"              "Papurana_jimiensis"          
[127] "Papurana_kreffti"             "Papurana_milleti"            
[129] "Papurana_milneana"            "Papurana_moluccana"          
[131] "Papurana_novaeguineae"        "Papurana_papua"              
[133] "Papurana_supragrisea"         "Papurana_volkerjane"         
[135] "Papurana_waliesa"             "Phrynoidis_asper"            
[137] "Pithecopus_hypochondrialis"   "Platymantis_banahao"         
[139] "Platymantis_bayani"           "Platymantis_biak"            
[141] "Platymantis_cagayanensis"     "Platymantis_cornutus"        
[143] "Platymantis_corrugatus"       "Platymantis_diesmosi"        
[145] "Platymantis_dorsalis"         "Platymantis_guentheri"       
[147] "Platymantis_hazelae"          "Platymantis_indeprensus"     
[149] "Platymantis_isarog"           "Platymantis_lawtoni"         
[151] "Platymantis_levigatus"        "Platymantis_luzonensis"      
[153] "Platymantis_mimulus"          "Platymantis_montanus"        
[155] "Platymantis_naomii"           "Platymantis_negrosensis"     
[157] "Platymantis_paengi"           "Platymantis_panayensis"      
[159] "Platymantis_polillensis"      "Platymantis_pseudodorsalis"  
[161] "Platymantis_pygmaeus"         "Platymantis_rabori"          
[163] "Platymantis_sierramadrensis"  "Platymantis_spelaeus"        
[165] "Platymantis_subterrestris"    "Platymantis_taylori"         
[167] "Plethodontohyla_alluaudi"     "Pulchrana_baramica"          
[169] "Pulchrana_debussyi"           "Pulchrana_glandulosa"        
[171] "Pulchrana_grandocula"         "Pulchrana_laterimaculata"    
[173] "Pulchrana_mangyanum"          "Pulchrana_melanomenta"       
[175] "Pulchrana_moellendorffi"      "Pulchrana_picturata"         
[177] "Pulchrana_siberu"             "Pulchrana_signata"           
[179] "Pulchrana_similis"            "Rentapia_everetti"           
[181] "Rentapia_hosii"               "Rhacophorus_rhyssocephalus"  
[183] "Rohanixalus_hansenae"         "Rohanixalus_vittatus"        
[185] "Sanguirana_mearnsi"           "Sclerophrys_buchneri"        
[187] "Sclerophrys_camerunensis"     "Sclerophrys_funerea"         
[189] "Sclerophrys_gracilipes"       "Sclerophrys_gutturalis"      
[191] "Sclerophrys_kassasii"         "Sclerophrys_latifrons"       
[193] "Sclerophrys_maculata"         "Sclerophrys_pentoni"         
[195] "Sclerophrys_regularis"        "Sclerophrys_superciliaris"   
[197] "Sclerophrys_tihamica"         "Sclerophrys_togoensis"       
[199] "Sclerophrys_tuberosa"         "Sclerophrys_xeros"           
[201] "Sphenophryne_allisoni"        "Sphenophryne_brevicrus"      
[203] "Sphenophryne_coggeri"         "Sphenophryne_crassa"         
[205] "Sphenophryne_dentata"         "Sphenophryne_magnitympanum"  
[207] "Sphenophryne_rhododactyla"    "Sphenophryne_rubra"          
[209] "Sphenophryne_schlaginhaufeni" "Sphenophryne_similis"        
[211] "Sphenophryne_stenodactyla"    "Sphenophryne_thomsoni"       
[213] "Sylvirana_cubitalis"          "Sylvirana_guentheri"         
[215] "Sylvirana_maosonensis"        "Sylvirana_mortenseni"        
[217] "Sylvirana_nigrovittata"       "Sylvirana_spinulosa"         
[219] "Trachycephalus_typhonius"     "Uperodon_nagaoi"             
[221] "Uperodon_obscurus"            "Uperodon_palmatus"           
[223] "Uperodon_taprobanicus"        "Uperodon_variegatus"         
[225] "Wijayarana_masonii"           "Wijayarana_melasma"          
[227] "Wijayarana_modiglianii"       "Wijayarana_sumatrana"        
[229] "Zhangixalus_achantharrhena"   "Zhangixalus_arboreus"        
[231] "Zhangixalus_arvalis"          "Zhangixalus_aurantiventris"  
[233] "Zhangixalus_dennysi"          "Zhangixalus_dulitensis"      
[235] "Zhangixalus_moltrechti"       "Zhangixalus_owstoni"         
[237] "Zhangixalus_prasinatus"       "Zhangixalus_prominanus"      
[239] "Zhangixalus_schlegelii"       "Zhangixalus_taipeianus"      
[241] "Zhangixalus_viridis" 