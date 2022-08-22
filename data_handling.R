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
#dir()
#anura_islands_orig <- as_tibble(read.csv2("comunidades_ilhas.csv", sep=","))
#head(anura_islands_orig)
#anura_islands <- anura_islands_orig %>% select(-X, -OBJECTID_1) #copy

# Sun Aug 14 11:00:01 2022 ------------------------------
# using dataset without 'ilhas fluviais'
dir()
matrix_islands <- as_tibble(read.csv2("first_dataset/matrix_islands_ok.csv", sep=",", fileEncoding='UTF-8'))
#removing name
matrix_islands$name_ID <- gsub("\\_..*","",matrix_islands$name_ID)

###
head(matrix_islands) # 2329
anura_islands_ok <- matrix_islands %>% select(-X, -name_ID)  

#check data
head(matrix_islands) 
head(anura_islands_ok) 

# we use copy data without island id columns 
nrow(anura_islands_ok) # 1499 islands 
ncol(anura_islands_ok) # 2327 anuran islands spp.

# traits data
dir()
anura_traits <- read.csv2("first_dataset/functional_anura.csv", sep=",") 
names(anura_traits) 
dim(anura_traits) #2249 spp and 34 traits

# selecting traits for our analysis
traits <- c("Species", "Ter","Fos","Aqu","Arb", "Dir", "Lar", "Viv", "Body_size_mm")
anura_traits_select <- anura_traits %>% select(traits)

# Thu Aug 11 12:05:59 2022 ------------------------------
# Join trait data with community data

#Conferindo quais espécies da lista possuem traits
# Lista para join
# Arrumando para linhas
View(anura_islands_ok)

anuran_list <- anura_islands_ok %>%
  gather(key="Species",value="count", Abavorana_luctuosa:Zhangixalus_yinggelingensis) %>% distinct(Species, .keep_all = TRUE)

nrow(anuran_list) #2327 spp composition matrix
nrow(anura_traits_select) # 2249 spp traits
# 78 spp losing in composition matrix

# Left join: espécies que estão no x mantidas e agrupa com os correspiondentes de y 
anura_with_traits <- left_join(anuran_list, anura_traits_select, by="Species") %>% select(-count)
vis_miss(anura_with_traits)
# habitat use (ter, fos, aqu, arb) multi-choice traits with na's
# breeding strategy (dir, lar, viv) choice trait - 10.89% missing data
# litter_size_max continuous trait - 68% missing data.  
# View(anura_traits_select)

# Anti join: espécies que estão no x, mas não estão no y (x,y)
anura_without_traits <- anti_join(anuran_list, anura_traits_select, by="Species")
nrow(anura_without_traits) # 247 spp.

#Retirando espécies da comunidade que não estão nos traits
rem.col.list <- anura_without_traits$Species
short_anura_islands <- anura_islands_ok[,!(names(anura_islands_ok)%in% rem.col.list)]
ncol(short_anura_islands) #2080 species with traits

###---
# phylogenetic data 
amphibia_phy <- read.tree("amph_shl_new_Consensus_7238.tre")

# we changed names in phylogeny for match to traits and composition matrix
name_phy <- as_tibble(amphibia_phy$tip.label) %>% dplyr::mutate(new_name = amphibia_phy$tip.label)
View(name_phy)

# nomenclature adjustments "Boana_semilineatus"
name_phy$new_name <- name_phy$new_name %>% str_replace_all("Hypsiboas", "Boana") %>% str_replace_all("Boana_albomarginatus", "Boana_albomarginata") %>% str_replace_all("Boana_calcaratus", "Boana_calcarata") %>% str_replace_all("Boana_fasciatus", "Boana_fasciata") %>% str_replace_all("Boana_geographicus", "Boana_geographica") %>% str_replace_all("Boana_pulchellus", "Boana_pulchella") %>% str_replace_all("Boana_multifasciatus", "Boana_multifasciata") %>% str_replace_all("Boana_ornatissimus", "Boana_ornatissima") %>% str_replace_all("Boana_punctatus", "Boana_punctata") %>% str_replace_all("Boana_rufitelus", "Boana_rufitela") %>% str_replace_all("Boana_semilineatus","Boana_semilineata") %>% str_replace_all("Boana_semilineatus", "Boana_semilineata") %>% str_replace_all("Chiromantis_doriae", "Chirixalus_doriae") %>% str_replace_all("Hylomantis_aspera", "Agalychnis_aspera") %>% str_replace_all("Pachymedusa_dacnicolor", "Agalychnis_dacnicolor") %>% str_replace_all("Hylarana_luctuosa", "Abavorana_luctuosa") %>% str_replace_all("Ingerana_baluensis", "Alcalus_baluensis") %>% str_replace_all("Ingerana_mariae", "Alcalus_mariae") %>% str_replace_all("Ingerana_rajae","Alcalus_rajae") %>% str_replace_all("Ingerana_sariba","Alcalus_sariba") %>% str_replace_all("Stumpffia_helenae","Anilany_helenae") %>% str_replace_all("Cophyla_phyllodactyla","Platypelis_phyllodactyla") %>% str_replace_all("Phyllomedusa_tomopterna", "Callimedusa_tomopterna") %>% str_replace_all("Syncope_hudsoni", "Chiasmocleis_hudsoni") %>% str_replace_all("Cyclorana_alboguttata", "Litoria_alboguttata") %>% str_replace_all("Cyclorana_australis", "Litoria_australis") %>% str_replace_all("Cyclorana_brevipes", "Litoria_brevipes")  %>% str_replace_all("Cyclorana_cryptotis", "Litoria_cryptotis") %>% str_replace_all("Cyclorana_longipes", "Litoria_longipes") %>% str_replace_all("Cyclorana_novaehollandiae", "Litoria_novaehollandiae") %>% str_replace_all("Cyclorana_nullicedens", "Litoria_nullicedens") %>% str_replace_all("Cyclorana_vagitus", "Litoria_vagitus") %>% str_replace_all("Scinax_trapicheiroi", "Ololygon_trapicheiroi") %>% str_replace_all("Scinax_agilis", "Ololygon_agilis") %>% str_replace_all("Scinax_argyreornatus", "Ololygon_argyreornata") %>% str_replace_all("Scinax_berthae", "Ololygon_berthae") %>% str_replace_all("Scinax_brieni", "Ololygon_brieni") %>% str_replace_all("Scinax_catharinae", "Ololygon_catharinae") %>% str_replace_all("Scinax_flavoguttatus", "Ololygon_flavoguttata") %>% str_replace_all("Scinax_humilis", "Ololygon_humilis") %>% str_replace_all("Scinax_jureia", "Ololygon_jureia") %>% str_replace_all("Scinax_littoralis", "Ololygon_littoralis") %>% str_replace_all("Scinax_perpusillus", "Ololygon_perpusilla") %>% str_replace_all("Scinax_rizibilis", "Ololygon_rizibilis") %>% str_replace_all("Rana_palmipes", "Lithobates_palmipes") %>% str_replace_all("Rana_berlandieri", "Lithobates_berlandieri") %>% str_replace_all("Rana_capito", "Lithobates_capito") %>% str_replace_all("Rana_catesbeiana", "Lithobates_catesbeianus") %>% str_replace_all("Rana_clamitans", "Lithobates_clamitans") %>% str_replace_all("Rana_forreri", "Lithobates_forreri") %>% str_replace_all("Rana_grylio", "Lithobates_grylio")  %>% str_replace_all("Rana_heckscheri", "Lithobates_heckscheri") %>% str_replace_all("Rana_magnaocularis", "Lithobates_magnaocularis")  %>% str_replace_all("Rana_palustris", "Lithobates_palustris")  %>% str_replace_all("Rana_pipiens", "Lithobates_pipiens")  %>% str_replace_all("Rana_septentrionalis", "Lithobates_septentrionalis") %>% str_replace_all("Rana_sphenocephala", "Lithobates_sphenocephalus") %>% str_replace_all("Rana_sylvatica", "Lithobates_sylvaticus") %>% str_replace_all("Rana_vaillanti", "Lithobates_vaillanti") %>% str_replace_all("Rana_virgatipes", "Lithobates_virgatipes")  %>% str_replace_all("Rana_warszewitschii", "Lithobates_warszewitschii") %>% str_replace_all("Rana_yavapaiensis", "Lithobates_yavapaiensis") %>% str_replace_all("Lysapsus_boliviana", "Lysapsus_bolivianus") 

name_phy$new_name <- name_phy$new_name %>% str_replace_all("Fejervarya_andamanensis", "Minervarya_andamanensis") %>% str_replace_all("Fejervarya_greenii", "Minervarya_greenii") %>% str_replace_all("Fejervarya_nepalensis", "Minervarya_nepalensis") %>% str_replace_all("Fejervarya_nicobariensis", "Minervarya_nicobariensis") %>% str_replace_all("Fejervarya_pierrei", "Minervarya_pierrei") %>% str_replace_all("Fejervarya_rufescens", "Minervarya_rufescens") %>% str_replace_all("Fejervarya_teraiensis", "Minervarya_teraiensis") %>% str_replace_all("Leptolalax_arayai", "Leptobrachella_arayai") %>% str_replace_all("Leptolalax_dringi", "Leptobrachella_dringi") %>% str_replace_all("Leptolalax_fuliginosus", "Leptobrachella_fuliginosa" ) %>% str_replace_all("Leptolalax_gracilis","Leptobrachella_gracilis") %>% str_replace_all("Leptolalax_hamidi","Leptobrachella_hamidi") %>% str_replace_all("Leptolalax_kajangensis","Leptobrachella_kajangensis") %>% str_replace_all("Leptolalax_liui","Leptobrachella_liui") %>% str_replace_all("Leptolalax_maurus", "Leptobrachella_maura") %>% str_replace_all("Leptolalax_pelodytoides","Leptobrachella_pelodytoides") %>% str_replace_all("Leptolalax_picta","Leptobrachella_pictus") %>% str_replace_all("Fejervarya_kirtisinghei", "Minervarya_kirtisinghei") %>% str_replace_all("Liophryne_allisoni","Sphenophryne_allisoni") %>% str_replace_all("Liophryne_dentata", "Sphenophryne_dentata") %>% str_replace_all("Liophryne_magnitympanum","Sphenophryne_magnitympanum") %>% str_replace_all("Liophryne_rhododactyla","Sphenophryne_rhododactyla") %>% str_replace_all("Liophryne_rubra","Sphenophryne_rubra") %>% str_replace_all("Liophryne_rubra","Sphenophryne_rubra") %>% str_replace_all("Liophryne_schlaginhaufeni", "Sphenophryne_schlaginhaufeni") %>% str_replace_all("Liophryne_similis","Sphenophryne_similis") %>% str_replace_all("Oxydactyla_brevicrus","Sphenophryne_brevicrus") %>% str_replace_all("Oxydactyla_coggeri","Sphenophryne_coggeri") %>% str_replace_all("Oxydactyla_crassa","Sphenophryne_crassa") %>% str_replace_all("Oxydactyla_stenodactyla","Sphenophryne_stenodactyla") %>% str_replace_all("Rhacophorus_achantharrhena", "Zhangixalus_achantharrhena") %>% str_replace_all("Rhacophorus_arboreus", "Zhangixalus_arboreus") %>% str_replace_all("Rhacophorus_arvalis", "Zhangixalus_arvalis") %>% str_replace_all("Rhacophorus_aurantiventris", "Zhangixalus_aurantiventris") %>% str_replace_all("Rhacophorus_dennysi", "Zhangixalus_dennysi") %>% str_replace_all("Rhacophorus_dulitensis", "Zhangixalus_dulitensis") %>% str_replace_all("Rhacophorus_moltrechti", "Zhangixalus_moltrechti") %>% str_replace_all("Rhacophorus_owstoni", "Zhangixalus_owstoni") %>% str_replace_all("Rhacophorus_prasinatus", "Zhangixalus_prasinatus") %>% str_replace_all("Rhacophorus_prominanus", "Zhangixalus_prominanus") %>% str_replace_all("Rhacophorus_schlegelii", "Zhangixalus_schlegelii") %>% str_replace_all("Rhacophorus_taipeianus", "Zhangixalus_taipeianus") %>% str_replace_all("Rhacophorus_viridis", "Zhangixalus_viridis") 

name_phy$new_name <- name_phy$new_name %>% str_replace_all("Hylarana_crassiovis", "Chalcorana_crassiovis") %>% str_replace_all("Hylarana_chalconota","Chalcorana_chalconota") %>% str_replace_all("Hylarana_kampeni","Chalcorana_kampeni") %>% str_replace_all("Hylarana_macrops","Chalcorana_macrops") %>% str_replace_all("Hylarana_megalonesa","Chalcorana_megalonesa") %>% str_replace_all("Hylarana_mocquardii", "Chalcorana_mocquardii") %>% str_replace_all("Hylarana_parvacola","Chalcorana_parvaccola") %>% str_replace_all("Hylarana_raniceps","Chalcorana_raniceps") %>% str_replace_all("Hylarana_rufipes", "Chalcorana_rufipes") %>% str_replace_all("Albericus_brunhildae","Choerophryne_brunhildae") %>% str_replace_all("Albericus_darlingtoni","Choerophryne_darlingtoni") %>% str_replace_all("Albericus_exclamitans","Choerophryne_exclamitans") %>% str_replace_all("Albericus_fafniri","Choerophryne_fafniri") %>% str_replace_all("Albericus_gunnari","Choerophryne_gunnari") %>% str_replace_all("Albericus_laurini","Choerophryne_laurini") %>% str_replace_all("Albericus_rhenaurum","Choerophryne_rhenaurum") %>% str_replace_all("Albericus_sanguinopictus","Choerophryne_sanguinopicta") %>% str_replace_all("Albericus_siegfriedi","Choerophryne_siegfriedi") %>% str_replace_all("Albericus_swanhildae","Choerophryne_swanhildae") %>% str_replace_all("Albericus_tuberculus","Choerophryne_tubercula") %>% str_replace_all("Albericus_valkuriarum","Choerophryne_valkuriarum") %>% str_replace_all("Albericus_variegatus","Choerophryne_variegata") %>% str_replace_all("Platypelis_barbouri","Cophyla_barbouri") %>% str_replace_all("Platypelis_grandis", "Cophyla_grandis") %>% str_replace_all("Platypelis_mavomavo", "Cophyla_mavomavo") %>% str_replace_all("Platypelis_milloti", "Cophyla_milloti") %>% str_replace_all("Platypelis_phyllodactyla", "Cophyla_phyllodactyla") %>% str_replace_all("Platypelis_pollicaris" ,"Cophyla_pollicaris") %>% str_replace_all("Platypelis_tetra", "Cophyla_tetra") %>% str_replace_all("Platypelis_tsaratananaensis", "Cophyla_tsaratananaensis") %>% str_replace_all("Platypelis_tuberifera", "Cophyla_tuberifera") %>% 
  str_replace_all("Platymantis_wuenscheorum","Cornufer_wuenscheorum") %>% str_replace_all("Batrachylodes_wolfi", "Cornufer_wolfi") %>% str_replace_all("Platymantis_weberi", "Cornufer_weberi") %>% str_replace_all("Discodeles_vogti", "Cornufer_vogti") %>% str_replace_all("Platymantis_vitiensis", "Cornufer_vitiensis") %>% str_replace_all("Platymantis_vitianus","Cornufer_vitianus") %>% str_replace_all("Batrachylodes_vertebralis","Cornufer_vertebralis") %>% str_replace_all("Batrachylodes_trossulus","Cornufer_trossulus") %>% str_replace_all("Platymantis_sulcatus","Cornufer_sulcatus") %>% str_replace_all("Platymantis_solomonis","Cornufer_solomonis") %>% str_replace_all("Platymantis_schmidti","Cornufer_schmidti") %>% str_replace_all("Platymantis_punctatus","Cornufer_punctatus") %>% str_replace_all("Platymantis_pelewensis","Cornufer_pelewensis") %>% str_replace_all("Platymantis_parkeri","Cornufer_parkeri") %>% str_replace_all("Platymantis_parilis","Cornufer_parilis") %>% str_replace_all("Platymantis_papuensis","Cornufer_papuensis") %>% str_replace_all("Discodeles_opisthodon","Cornufer_opisthodon") %>% str_replace_all("Platymantis_nexipus","Cornufer_nexipus") %>% str_replace_all("Platymantis_neckeri","Cornufer_neckeri") %>% str_replace_all("Platymantis_nakanaiorum","Cornufer_nakanaiorum") %>% str_replace_all("Platymantis_myersi","Cornufer_myersi") %>% str_replace_all("Batrachylodes_montanus","Cornufer_montanus") %>% str_replace_all("Batrachylodes_minutus", "Cornufer_minutus") %>% str_replace_all("Platymantis_mimicus","Cornufer_mimicus") %>% str_replace_all("Batrachylodes_mediodiscus","Cornufer_mediodiscus") %>% str_replace_all("Platymantis_manus","Cornufer_manus") %>% str_replace_all("Platymantis_mamusiorum","Cornufer_mamusiorum") %>% str_replace_all("Discodeles_malukuna","Cornufer_malukuna") %>% str_replace_all("Platymantis_magnus","Cornufer_magnus") %>% str_replace_all("Platymantis_macrosceles","Cornufer_macrosceles") 

name_phy$new_name <- name_phy$new_name %>% str_replace_all("Platymantis_macrops","Cornufer_macrops") %>% str_replace_all("Platymantis_latro","Cornufer_latro") %>% str_replace_all("Discodeles_guppyi","Cornufer_guppyi") %>% str_replace_all("Platymantis_guentheri","Cornufer_guentheri") %>% str_replace_all("Platymantis_gilliardi","Cornufer_gilliardi") %>% str_replace_all("Batrachylodes_gigas","Cornufer_gigas") %>% str_replace_all("Batrachylodes_elegans","Cornufer_elegans") %>% str_replace_all("Platymantis_desticans","Cornufer_desticans") %>% str_replace_all("Platymantis_cryptotis","Cornufer_cryptotis") %>% str_replace_all("Platymantis_cheesmanae","Cornufer_cheesmanae") %>% str_replace_all("Platymantis_caesiops","Cornufer_caesiops") %>% str_replace_all("Platymantis_bufonulus","Cornufer_bufonulus") %>% str_replace_all("Platymantis_browni","Cornufer_browni") %>% str_replace_all("Platymantis_boulengeri","Cornufer_boulengeri") %>% str_replace_all("Platymantis_bimaculatus","Cornufer_bimaculatus") %>% str_replace_all("Platymantis_batantae","Cornufer_batantae") %>% str_replace_all("Platymantis_akarithymus","Cornufer_akarithymus") %>% str_replace_all("Platymantis_admiraltiensis", "Cornufer_admiraltiensis") %>% str_replace_all("Platymantis_adiastolus","Cornufer_adiastolus") %>% str_replace_all("Platymantis_aculeodactylus", "Cornufer_aculeodactylus") %>% str_replace_all("Platymantis_acrochordus","Cornufer_acrochordus") %>% str_replace_all("Hyla_andersonii","Dryophytes_andersonii") %>% str_replace_all("Hyla_chrysoscelis","Dryophytes_chrysoscelis") %>% str_replace_all("Hyla_avivoca","Dryophytes_avivoca") %>% str_replace_all("Hyla_cinerea","Dryophytes_cinereus") %>% str_replace_all("Hyla_femoralis","Dryophytes_femoralis") %>% str_replace_all("Hyla_gratiosa","Dryophytes_gratiosus") %>% str_replace_all("Hyla_femoralis","Dryophytes_femoralis") %>% str_replace_all("Hyla_immaculata","Dryophytes_immaculatus") %>% str_replace_all("Hyla_japonica","Dryophytes_japonicus") %>% str_replace_all("Hyla_squirella","Dryophytes_squirellus") %>% str_replace_all("Hyla_versicolor","Dryophytes_versicolor") %>% str_replace_all("Borneophrys_edwardinae", "Megophrys_edwardinae") %>% str_replace_all("Xenophrys_aceras", "Megophrys_aceras") %>% str_replace_all("Xenophrys_baluensis", "Megophrys_baluensis") %>% str_replace_all("Xenophrys_boettgeri", "Megophrys_boettgeri") %>% str_replace_all("Xenophrys_brachykolos", "Megophrys_brachykolos") %>% str_replace_all("Xenophrys_dringi", "Megophrys_dringi") %>% str_replace_all("Xenophrys_major", "Megophrys_major") %>% str_replace_all("Xenophrys_parallela", "Megophrys_parallela") %>%  str_replace_all("Xenophrys_parva", "Megophrys_parva") %>% str_replace_all("Hylarana_attigua","Papurana_attigua") %>% str_replace_all("Hylarana_aurata","Papurana_aurata") %>% str_replace_all("Hylarana_daemeli","Papurana_daemeli") %>% str_replace_all("Hylarana_elberti","Papurana_elberti") %>% str_replace_all("Hylarana_florensis","Papurana_florensis") %>% str_replace_all("Hylarana_garritor","Papurana_garritor") %>% str_replace_all("Hylarana_grisea","Papurana_grisea") %>% str_replace_all("Hylarana_jimiensis","Papurana_jimiensis") %>% str_replace_all("Hylarana_kreffti","Papurana_kreffti") %>% str_replace_all("Hylarana_milleti","Papurana_milleti") %>% str_replace_all("Hylarana_milneana","Papurana_milneana") %>% str_replace_all("Hylarana_moluccana","Papurana_moluccana") %>% str_replace_all("Hylarana_novaeguineae","Papurana_novaeguineae") %>% str_replace_all("Hylarana_papua","Papurana_papua") %>% str_replace_all("Hylarana_supragrisea","Papurana_supragrisea") %>% str_replace_all("Hylarana_volkerjane","Papurana_volkerjane") %>% str_replace_all("Hylarana_waliesa","Papurana_waliesa") 

name_phy$new_name <- name_phy$new_name %>% str_replace_all("Amietophrynus_buchneri","Sclerophrys_buchneri") %>% str_replace_all("Amietophrynus_camerunensis", "Sclerophrys_camerunensis") %>% str_replace_all("Amietophrynus_funereus","Sclerophrys_funerea") %>% str_replace_all("Amietophrynus_gracilipes","Sclerophrys_gracilipes") %>% str_replace_all("Amietophrynus_gutturalis","Sclerophrys_gutturalis") %>% str_replace_all("Amietophrynus_kassasii","Sclerophrys_kassasii") %>% str_replace_all("Amietophrynus_latifrons","Sclerophrys_latifrons") %>% str_replace_all("Amietophrynus_maculatus","Sclerophrys_maculata") %>% str_replace_all("Amietophrynus_regularis","Sclerophrys_regularis") %>% str_replace_all("Amietophrynus_superciliaris","Sclerophrys_superciliaris") %>% str_replace_all("Amietophrynus_togoensis","Sclerophrys_togoensis" ) %>% str_replace_all("Amietophrynus_tuberosus","Sclerophrys_tuberosa") %>% str_replace_all("Amietophrynus_xeros","Sclerophrys_xeros") %>% str_replace_all("Calluella", "Glyphoglossus") %>% str_replace_all("Glyphoglossus_guttulata", "Glyphoglossus_guttulatus") %>% str_replace_all("Glyphoglossus_flava", "Glyphoglossus_flavus") %>% str_replace_all("Hylarana_albolabris","Amnirana_albolabris") %>% str_replace_all("Hylarana_amnicola","Amnirana_amnicola") %>% str_replace_all("Hylarana_galamensis","Amnirana_galamensis") %>% str_replace_all("Hylarana_lepus","Amnirana_lepus") %>% str_replace_all("Hylarana_occidentalis","Amnirana_occidentalis") %>% str_replace_all("Hylarana_baramica","Pulchrana_baramica") %>% str_replace_all("Hylarana_debussyi","Pulchrana_debussyi") %>% str_replace_all("Hylarana_glandulosa","Pulchrana_glandulosa") %>% str_replace_all("Hylarana_grandocula","Pulchrana_grandocula") %>% str_replace_all("Hylarana_laterimaculata","Pulchrana_laterimaculata") %>% str_replace_all("Hylarana_mangyanum","Pulchrana_mangyanum") %>% str_replace_all("Hylarana_melanomenta","Pulchrana_melanomenta") %>% str_replace_all("Hylarana_moellendorffi","Pulchrana_moellendorffi") %>% str_replace_all("Hylarana_picturata","Pulchrana_picturata") %>% str_replace_all("Hylarana_siberu","Pulchrana_siberu") %>% str_replace_all("Hylarana_signata","Pulchrana_signata") %>% str_replace_all("Hylarana_similis","Pulchrana_similis")%>% str_replace_all("Oreophryne_nana","Aphantophryne_nana") %>% str_replace_all("Pseudocallulops_eurydactylus","Asterophrys_eurydactyla") %>% str_replace_all("Pseudocallulops_pullifer","Asterophrys_pullifer") %>% str_replace_all("Metamagnusia_marani","Asterophrys_marani") %>% str_replace_all("Metamagnusia_slateri","Asterophrys_slateri") %>% str_replace_all("Hylarana_nicobariensis","Bijurana_nicobariensis") %>% str_replace_all("Chiromantis_nongkhorensis","Chirixalus_nongkhorensis") %>% str_replace_all("Oxydactyla_alpestris","Copiula_alpestris") %>% str_replace_all("Austrochaperina_guttata","Copiula_guttata") %>% str_replace_all("Huia_masonii","Wijayarana_masonii") %>% str_replace_all("Huia_melasma","Wijayarana_melasma") %>% str_replace_all("Huia_modiglianii","Wijayarana_modiglianii") %>% str_replace_all("Huia_sumatrana","Wijayarana_sumatrana") %>% str_replace_all("Rugosa_emelianjovi","Glandirana_emeljanovi") %>% str_replace_all("Rugosa_rugosa","Glandirana_rugosa") %>% str_replace_all("Ramanella_obscura","Uperodon_obscurus") %>% str_replace_all("	
Ramanella_nagaoi","Uperodon_nagaoi") %>% str_replace_all("Ramanella_palmata","Uperodon_palmatus") %>% str_replace_all("Ramanella_variegata","Uperodon_variegatus") %>% str_replace_all("Kaloula_taprobanica","Uperodon_taprobanicus") %>% str_replace_all("Hylarana_cubitalis","Sylvirana_cubitalis") %>% str_replace_all("Hylarana_guentheri","Sylvirana_guentheri") %>% str_replace_all("Hylarana_maosonensis","Sylvirana_maosonensis") %>% str_replace_all("Hylarana_mortenseni","Sylvirana_mortenseni") %>% str_replace_all("Hylarana_nigrovittata","Sylvirana_nigrovittata") %>% str_replace_all("Hylarana_spinulosa","Sylvirana_spinulosa") %>% str_replace_all("Phyllomedusa_hypochondrialis","Pithecopus_hypochondrialis") 

name_phy$new_name <- name_phy$new_name %>% str_replace_all("Hylarana_malabarica","Hydrophylax_malabaricus") %>% str_replace_all("Hylarana_leptoglossa","Hydrophylax_leptoglossa") %>% str_replace_all("Hylarana_gracilis","Hydrophylax_gracilis") %>% str_replace_all("Mantophryne_infulata","Hylophorbus_infulatus") %>% str_replace_all("Chlorolius_koehleri","Hyperolius_koehleri") %>% str_replace_all("Hylarana_aurantiaca","Indosylvirana_aurantiaca") %>% str_replace_all("Hylarana_temporalis","Indosylvirana_temporalis") %>% str_replace_all("Leptolalax_pictus","Leptobrachella_picta") %>% str_replace_all("Limnonectes_rhacoda","Limnonectes_rhacodus") %>% str_replace_all("Nyctimystes_rueppelli","Litoria_rueppelli") %>% str_replace_all("Pherohapsis_menziesi","Mantophryne_menziesi") %>% str_replace_all("Microhyla_erythropoda","Micryletta_erythropoda") %>% str_replace_all("Fejervarya_syhadrensis","Minervarya_syhadrensis") %>% str_replace_all("Microhyla_perparva","Nanohyla_perparva") %>% str_replace_all("Microhyla_petrigena","Nanohyla_petrigena") %>% str_replace_all("Babina_adenopleura","Nidirana_adenopleura") %>% str_replace_all("Babina_hainanensis","Nidirana_hainanensis") %>% str_replace_all("Babina_okinavana","Nidirana_okinavana") %>% str_replace_all("Litoria_infrafrenata","Nyctimystes_infrafrenatus") %>% str_replace_all("Litoria_sauroni","Nyctimystes_sauroni") %>% str_replace_all("Nyctimystes_tyleri","Litoria_tyleri") %>% str_replace_all("Hylarana_arfaki","Papurana_arfaki") %>% str_replace_all("Ramanella_nagaoi","Uperodon_nagaoi") %>% str_replace_all("Phrynoidis_aspera","Phrynoidis_asper") %>% str_replace_all("Cornufer_guentheri","Platymantis_guentheri") %>% str_replace_all("Rhombophryne_alluaudi","Plethodontohyla_alluaudi") %>% str_replace_all("Pedostibes_everetti","Rentapia_everetti") %>% str_replace_all("Pedostibes_hosii","Rentapia_hosii") %>% str_replace_all("Chiromantis_hansenae","Rohanixalus_hansenae") %>% str_replace_all("Chiromantis_vittatus","Rohanixalus_vittatus") %>% str_replace_all("Genyophryne_thomsoni","Sphenophryne_thomsoni")

amphibia_phy$tip.label <- name_phy$new_name
View(name_phy)
# Checking species match
options(max.print=200)
match.phylo.comm(amphibia_phy,short_anura_islands)

#Prune with phylogeny Pyron and Wiens consensus
phy_anura_islands <- prune.sample(short_anura_islands, amphibia_phy) 

phy_anura_islands # 2068 spp

ncol(short_anura_islands) #2080
#quantas espécies perdemos?
2080-2068 # 12 spp

# Phylogeny anuran islands
plot(phy_anura_islands, type="fan",show.tip.label = FALSE)
#save phy islands
write.tree(phy_anura_islands, file = "anura_islands.tre")

###
# "Anodonthyla_boulengeri","Megophrys_carinense" nomes inválidos segundo Frost

# espécie duplicada na composição e nos traits  e litoria tyleri
missing.in.phy <- c("Anodonthyla_boulengeri","Cornufer_guentheri","Cornufer_hedigeri","Cornufer_heffernani","Fejervarya_pulla","Megophrys_carinense","Nyctimystes_tyleri","Rhacophorus_rhyssocephalus","Sanguirana_mearnsi", "Sclerophrys_pentoni","Sclerophrys_tihamica","Trachycephalus_typhonius")

phy_community_anura <- short_anura_islands[,!(names(short_anura_islands)%in% missing.in.phy)]

# Join lista da comunidade e filogenia com os traits
# Arrumando para linhas
anuran_list_phy <- phy_community_anura %>%
  gather(key="Species",value="count", Abavorana_luctuosa:Zhangixalus_viridis) %>% distinct(Species, .keep_all = TRUE)
nrow(anuran_list_phy) #2068 spp composição de espécies
nrow(anura_traits_select) # 2249 spp traits

# Left join: espécies que estão no x mantidas e agrupa com os correspiondentes de y 
phy_community_traits <- left_join(anuran_list_phy, anura_traits_select, by="Species") %>% select(-count)
vis_miss(anura_with_traits)

#Data
ncol(phy_community_anura) #2068 species community
phy_anura_islands #2068 species phylogeny
nrow(phy_community_traits) #2068 species traits

# Mon Aug 22 13:47:24 2022 ------------------------------
# taking off offspring size and adding body size

# traits raoni
traits_anura_raoni <- openxlsx::read.xlsx("first_dataset/functional_traits_raoni.xlsx")
names(traits_anura_raoni)
traits_anura_raoni <- as_tibble(traits_anura_raoni) %>% select("Species","Body_size_mm", "Fos","Ter","Aqu","Arb", "Dir","Lar","Viv")

# Join with list data
traits_anura_islands <- left_join(anuran_list_phy, traits_anura_raoni, by="Species") %>% select(-count)
vis_miss(traits_anura_islands)

# Saving traits data
write.csv2(phy_community_traits, "anura_islands_traits.csv", sep=",")
write.csv2(phy_community_anura, "anura_islands_matrix.csv", sep=",")
View(phy_community_anura)

# Mon Aug 22 14:00:28 2022 ------------------------------
# saving traits_anura_islands
dir()
write.csv2(traits_anura_islands, "anura_traits_raoni.csv", sep=",")

# Sun Aug 14 16:33:06 2022 ------------------------------
community_with_id <- mutate(phy_community_anura,id = matrix_islands$name_ID,)

# islands without anurans
index <- rowSums(phy_community_anura) #checking islands without anura spp
View(as.data.frame(index))

rem.lines <- c("17651", "18539", "18930", "22423", "3381", "3663", "5167", "9445")

community_without_zero <- community_with_id[!(community_with_id$id %in% rem.lines), ]

#Save
View(community_without_zero)
write.csv2(community_without_zero, "anura_islands_without0.csv", sep=",")