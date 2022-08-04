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
ncol(anura_islands) # 2503 anuran islands spp.

# phylogenetic data 
amphibia_phy <- read.tree("amph_shl_new_Consensus_7238.tre")
name_phy <- as.data.frame(amphibia_phy$tip.label)
View(name_phy)

# Checking species match
match.phylo.comm(amphibia_phy,anura_islands)
options(max.print=315)

# nomenclature adjustments 
# genus
names(anura_islands) <- gsub("Hylomantis","Agalychnis", fixed=T, names(anura_islands)) 
names(anura_islands) <- gsub("Boana", "Hypsiboas", fixed=T, names(anura_islands)) 
names(anura_islands) <- gsub("Amnirana", "Hylarana", fixed=T, names(anura_islands)) 
names(anura_islands) <- gsub("Cophyla", "Platypelis", fixed=T, names(anura_islands)) 

# epitete
names(anura_islands) <- gsub("Hypsiboas_albomarginata", "Hypsiboas_albomarginatus", fixed=T, names(anura_islands)) 
names(anura_islands) <- gsub("Hypsiboas_calcarata", "Hypsiboas_calcaratus", fixed=T, names(anura_islands)) 
names(anura_islands) <- gsub("Hypsiboas_fasciata", "Hypsiboas_fasciatus", fixed=T, names(anura_islands)) 
names(anura_islands) <- gsub("Hypsiboas_geographica", "Hypsiboas_geographicus", fixed=T, names(anura_islands)) 
names(anura_islands) <- gsub("Hypsiboas_multifasciata", "Hypsiboas_multifasciatus", fixed=T, names(anura_islands)) 
names(anura_islands) <- gsub("Hypsiboas_ornatissima", "Hypsiboas_ornatissimus", fixed=T, names(anura_islands)) 
names(anura_islands) <- gsub("Hypsiboas_punctata", "Hypsiboas_punctatus", fixed=T, names(anura_islands)) 
names(anura_islands) <- gsub("Hypsiboas_rufitela", "Hypsiboas_rufitelus", fixed=T, names(anura_islands)) 
names(anura_islands) <- gsub("Hypsiboas_semilineata", "Hypsiboas_semilineatus", fixed=T, names(anura_islands)) 
names(anura_islands) <- gsub("Chirixalus_doriae", "Chiromantis_doriae", fixed=T, names(anura_islands)) 
names(anura_islands) <- gsub("Agalychnis_aspera", "Hylomantis_aspera", fixed=T, names(anura_islands))
names(anura_islands) <- gsub("Agalychnis_dacnicolor", "Pachymedusa_dacnicolor", fixed=T, names(anura_islands))
names(anura_islands) <- gsub("Abavorana_luctuosa", "Hylarana_luctuosa", fixed=T, names(anura_islands))
names(anura_islands) <- gsub("Alcalus_baluensis" , "Ingerana_baluensis", fixed=T, names(anura_islands))
names(anura_islands) <- gsub("Alcalus_mariae"  , "Ingerana_mariae", fixed=T, names(anura_islands))
names(anura_islands) <- gsub("Alcalus_rajae"  ,"Ingerana_rajae", fixed=T, names(anura_islands))
names(anura_islands) <- gsub("Alcalus_sariba"  ,"Ingerana_sariba", fixed=T, names(anura_islands))
names(anura_islands) <- gsub("Anilany_helenae","Stumpffia_helenae", fixed=T, names(anura_islands))
names(anura_islands) <- gsub("Platypelis_phyllodactyla","Cophyla_phyllodactyla", fixed=T, names(anura_islands)) 

# Prune with phylogeny Pyron and Wiens consensus
phy_anura <- prune.sample(anura_islands, amphibia_phy) # 2051 spp

phy_anura 
anura_islands

# spp
[1] "Dropping taxa from the community because they are not present in the phylogeny:"                 
[3] "Aglyptodactylus_australis"     "Aglyptodactylus_chorus"       
[5] "Aglyptodactylus_inguinalis"    "Amolops_albispinus"           
[7] "Anodonthyla_boulengeri"        "Ansonia_phuketensis"          
[9] "Ansonia_teneritas"             "Ansonia_vidua"                
[11] "Aphantophryne_nana"            "Asterophrys_eurydactyla"      
[13] "Asterophrys_foja"              "Asterophrys_marani"           
[15] "Asterophrys_pullifer"          "Asterophrys_slateri"          
[17] "Austrochaperina_alexanderi"    "Austrochaperina_laurae"       
[19] "Austrochaperina_punctata"      "Austrochaperina_rudolfarndti" 
[21] "Bijurana_nicobariensis"        "Blommersia_nataliae"          
[23] "Blommersia_transmarina"        "Hypsiboas_greographicus"      
[25] "Hypsiboas_pulchella"           "Boophis_ankarafensis"         
[27] "Boophis_boppa"                 "Boophis_nauticus"             
[29] "Buergeria_choui"               "Buergeria_otai"               
[31] "Callimedusa_tomopterna"        "Callulops_argus"              
[33] "Callulops_biakensis"           "Callulops_bicolor"            
[35] "Callulops_neuhaussi"           "Callulops_stellatus"          
[37] "Callulops_taxispilotus"        "Callulops_wondiwoiensis"      
[39] "Callulops_yapenensis"          "Chalcorana_chalconota"        
[41] "Chalcorana_crassiovis"         "Chalcorana_kampeni"           
[43] "Chalcorana_macrops"            "Chalcorana_megalonesa"        
[45] "Chalcorana_mocquardii"         "Chalcorana_parvaccola"        
[47] "Chalcorana_raniceps"           "Chalcorana_rufipes"           
[49] "Chiasmocleis_hudsoni"          "Chirixalus_nongkhorensis"     
[51] "Chirixalus_trilaksonoi"        "Choerophryne_alainduboisi"    
[53] "Choerophryne_alpestris"        "Choerophryne_bickfordi"       
[55] "Choerophryne_bisyllaba"        "Choerophryne_brevicrus"       
[57] "Choerophryne_brunhildae"       "Choerophryne_crucifer"        
[59] "Choerophryne_darlingtoni"      "Choerophryne_epirrhina"       
[61] "Choerophryne_exclamitans"      "Choerophryne_fafniri"         
[63] "Choerophryne_grylloides"       "Choerophryne_gudrunae"        
[65] "Choerophryne_gunnari"          "Choerophryne_laurini"         
[67] "Choerophryne_multisyllaba"     "Choerophryne_murrita"         
[69] "Choerophryne_pandanicola"      "Choerophryne_pipiens"         
[71] "Choerophryne_rhenaurum"        "Choerophryne_sanguinopicta"   
[73] "Choerophryne_siegfriedi"       "Choerophryne_swanhildae"      
[75] "Choerophryne_tubercula"        "Choerophryne_valkuriarum"     
[77] "Choerophryne_variegata"        "Cophixalus_cateae"            
[79] "Cophixalus_hannahae"           "Cophixalus_rajampatensis"     
[81] "Cophixalus_salawatiensis"      "Cophixalus_viridis"           
[83] "Platypelis_berara"             "Platypelis_karenae"           
[85] "Platypelis_maharipeo"          "Platypelis_noromalalae"       
[87] "Platypelis_occultans"          "Platypelis_puellarum"         
[89] "Platypelis_rava"               "Copiula_alpestris"            
[91] "Copiula_annanoreenae"          "Copiula_guttata"              
[93] "Copiula_lennarti"              "Cornufer_acrochordus"         
[95] "Cornufer_aculeodactylus"       "Cornufer_adiastolus"          
[97] "Cornufer_admiraltiensis"       "Cornufer_akarithymus"         
[99] "Cornufer_batantae"             "Cornufer_bimaculatus"         
[101] "Cornufer_boulengeri"           "Cornufer_browni"              
[103] "Cornufer_bufoniformis"         "Cornufer_bufonulus"           
[105] "Cornufer_caesiops"             "Cornufer_cheesmanae"          
[107] "Cornufer_citrinospilus"        "Cornufer_cryptotis"           
[109] "Cornufer_custos"               "Cornufer_desticans"           
[111] "Cornufer_elegans"              "Cornufer_exedrus"             
[113] "Cornufer_gigas"                "Cornufer_gilliardi"           
[115] "Cornufer_guentheri"            "Cornufer_guppyi"              
[117] "Cornufer_hedigeri"             "Cornufer_heffernani"          
[119] "Cornufer_latro"                "Cornufer_macrops"             
[121] "Cornufer_macrosceles"          "Cornufer_magnus"              
[123] "Cornufer_malukuna"             "Cornufer_mamusiorum"          
[125] "Cornufer_manus"                "Cornufer_mediodiscus"         
[127] "Cornufer_mimicus"              "Cornufer_minutus"             
[129] "Cornufer_montanus"             "Cornufer_myersi"              
[131] "Cornufer_nakanaiorum"          "Cornufer_neckeri"             
[133] "Cornufer_nexipus"              "Cornufer_opisthodon"          
[135] "Cornufer_papuensis"            "Cornufer_parilis"             
[137] "Cornufer_parkeri"              "Cornufer_pelewensis"          
[139] "Cornufer_punctatus"            "Cornufer_schmidti"            
[141] "Cornufer_solomonis"            "Cornufer_sulcatus"            
[143] "Cornufer_trossulus"            "Cornufer_vertebralis"         
[145] "Cornufer_vitianus"             "Cornufer_vitiensis"           
[147] "Cornufer_vogti"                "Cornufer_weberi"              
[149] "Cornufer_wolfi"                "Cornufer_wuenscheorum"        
[151] "Dendropsophus_baileyi"         "Dryophytes_andersonii"        
[153] "Dryophytes_avivoca"            "Dryophytes_chrysoscelis"      
[155] "Dryophytes_cinereus"           "Dryophytes_femoralis"         
[157] "Dryophytes_gratiosus"          "Dryophytes_immaculatus"       
[159] "Dryophytes_japonicus"          "Dryophytes_squirellus"        
[161] "Dryophytes_versicolor"         "Eleutherodactylus_campi"      
[163] "Fejervarya_pulla"              "Glandirana_emeljanovi"        
[165] "Glandirana_rugosa"             "Glandirana_susurra"           
[167] "Glyphoglossus_brooksii"        "Glyphoglossus_capsus"         
[169] "Glyphoglossus_flavus"          "Glyphoglossus_guttulatus"     
[171] "Glyphoglossus_smithi"          "Glyphoglossus_volzi"          
[173] "Guibemantis_diphonus"          "Hydrophylax_gracilis"         
[175] "Hydrophylax_leptoglossa"       "Hydrophylax_malabaricus"      
[177] "Hyla_perrini"                  "Hylophorbus_infulatus"        
[179] "Hylophorbus_sigridae"          "Hyperolius_koehleri"          
[181] "Hyperolius_stictus"            "Hypopachus_ustus"             
[183] "Indosylvirana_aurantiaca"      "Indosylvirana_serendipi"      
[185] "Indosylvirana_temporalis"      "Kalophrynus_barioensis"       
[187] "Lankanectes_pera"              "Leptobrachella_arayai"        
[189] "Leptobrachella_dringi"         "Leptobrachella_fritinniens"   
[191] "Leptobrachella_fuliginosa"     "Leptobrachella_gracilis"      
[193] "Leptobrachella_hamidi"         "Leptobrachella_itiokai"       
[195] "Leptobrachella_juliandringi"   "Leptobrachella_kajangensis"   
[197] "Leptobrachella_laui"           "Leptobrachella_liui"          
[199] "Leptobrachella_marmorata"      "Leptobrachella_maura"         
[201] "Leptobrachella_pelodytoides"   "Leptobrachella_picta"         
[203] "Leptobrachella_sabahmontana"   "Leptobrachium_kantonishikawai"
[205] "Leptodactylus_pallidirostris"  "Leptophryne_javanica"         
[207] "Limnonectes_cintalubang"       "Limnonectes_larvaepartus"     
[209] "Limnonectes_quangninhensis"    "Limnonectes_rhacodus"         
[211] "Lithobates_berlandieri"        "Lithobates_capito"            
[213] "Lithobates_catesbeianus"       "Lithobates_clamitans"         
[215] "Lithobates_forreri"            "Lithobates_grylio"            
[217] "Lithobates_heckscheri"         "Lithobates_magnaocularis"     
[219] "Lithobates_palmipes"           "Lithobates_palustris"         
[221] "Lithobates_pipiens"            "Lithobates_septentrionalis"   
[223] "Lithobates_sphenocephalus"     "Lithobates_sylvaticus"        
[225] "Lithobates_vaillanti"          "Lithobates_virgatipes"        
[227] "Lithobates_warszewitschii"     "Lithobates_yavapaiensis"      
[229] "Litoria_alboguttata"           "Litoria_australis"            
[231] "Litoria_brevipes"              "Litoria_cryptotis"            
[233] "Litoria_longipes"              "Litoria_novaehollandiae"      
[235] "Litoria_nullicedens"           "Litoria_pallidofemora"        
[237] "Litoria_pinocchio"             "Litoria_rueppelli"            
[239] "Litoria_vagitus"               "Liuixalus_calcarius"          
[241] "Lysapsus_bolivianus"           "Mantophryne_insignis"         
[243] "Mantophryne_menziesi"          "Megophrys_aceras"             
[245] "Megophrys_baluensis"           "Megophrys_boettgeri"          
[247] "Megophrys_brachykolos"         "Megophrys_carinense"          
[249] "Megophrys_dringi"              "Megophrys_edwardinae"         
[251] "Megophrys_major"               "Megophrys_microstoma"         
[253] "Megophrys_parallela"           "Megophrys_parva"              
[255] "Meristogenys_penrissenensis"   "Microhyla_kuramotoi"          
[257] "Microhyla_mihintalei"          "Micryletta_erythropoda"       
[259] "Micryletta_nigromaculata"      "Minervarya_andamanensis"      
[261] "Minervarya_chiangmaiensis"     "Minervarya_greenii"           
[263] "Minervarya_kirtisinghei"       "Minervarya_muangkanensis"     
[265] "Minervarya_nepalensis"         "Minervarya_nicobariensis"     
[267] "Minervarya_pierrei"            "Minervarya_rufescens"         
[269] "Minervarya_syhadrensis"        "Minervarya_teraiensis"        
[271] "Nanohyla_perparva"             "Nanohyla_petrigena"           
[273] "Nidirana_adenopleura"          "Nidirana_hainanensis"         
[275] "Nidirana_okinavana"            "Nyctimystes_infrafrenatus"    
[277] "Nyctimystes_myolae"            "Nyctimystes_ocreptus"         
[279] "Nyctimystes_sauroni"           "Nyctimystes_tyleri"           
[281] "Ololygon_agilis"               "Ololygon_argyreornata"        
[283] "Ololygon_berthae"              "Ololygon_brieni"              
[285] "Ololygon_catharinae"           "Ololygon_flavoguttata"        
[287] "Ololygon_humilis"              "Ololygon_jureia"              
[289] "Ololygon_littoralis"           "Ololygon_perpusilla"          
[291] "Ololygon_rizibilis"            "Ololygon_trapicheiroi"        
[293] "Oreophryne_albitympanum"       "Oreophryne_albomaculata"      
[295] "Oreophryne_anser"              "Oreophryne_aurora"            
[297] "Oreophryne_banshee"            "Oreophryne_brunnea"           
[299] "Oreophryne_cameroni"           "Oreophryne_choerophrynoides"
[301] "Oreophryne_curator"            "Oreophryne_equus"             
[303] "Oreophryne_flavomaculata"      "Oreophryne_gagneorum"         
[305] "Oreophryne_graminis"           "Oreophryne_lemur"             
[307] "Oreophryne_matawan"            "Oreophryne_meliades"          
[309] "Oreophryne_nicolasi"           "Oreophryne_parkopanorum"      
[311] "Oreophryne_penelopeia"         "Oreophryne_philosylleptoris"  
[313] "Oreophryne_phoebe"             "Oreophryne_picticrus"         
[315] "Oreophryne_pseudunicolor"         "Oreophryne_roedeli"                
[317] "Paedophryne_titan"                  "Papurana_arfaki"                   
[319] "Papurana_attigua"                   "Papurana_aurata"                   
[321] "Papurana_daemeli"                   "Papurana_elberti"                  
[323] "Papurana_florensis"                 "Papurana_garritor"                 
[325] "Papurana_grisea"                    "Papurana_jimiensis"                
[327] "Papurana_kreffti"                   "Papurana_milleti"                  
[329] "Papurana_milneana"                  "Papurana_moluccana"                
[331] "Papurana_novaeguineae"              "Papurana_papua"                    
[333] "Papurana_supragrisea"               "Papurana_volkerjane"               
[335] "Papurana_waliesa"                   "Pelobates_balcanicus"              
[337] "Pelophryne_penrissenensis"          "Petropedetes_newtonii"             
[339] "Philautus_catbaensis"               "Philautus_nephophilus"             
[341] "Phrynoidis_asper"                   "Pithecopus_hypochondrialis"        
[343] "Platymantis_quezoni"                "Plethodontohyla_alluaudi"          
[345] "Polypedates_pseudotilophus"         "Pseudophilautus_conniffae"         
[347] "Pulchrana_baramica"                 "Pulchrana_centropeninsularis"      
[349] "Pulchrana_debussyi"                 "Pulchrana_glandulosa"              
[351] "Pulchrana_grandocula"               "Pulchrana_guttmani"                
[353] "Pulchrana_laterimaculata"           "Pulchrana_mangyanum"               
[355] "Pulchrana_melanomenta"              "Pulchrana_moellendorffi"           
[357] "Pulchrana_picturata"                "Pulchrana_rawa"                    
[359] "Pulchrana_siberu"                   "Pulchrana_signata"                 
[361] "Pulchrana_similis"                  "Rana_neba"                         
[363] "Rentapia_everetti"                  "Rentapia_hosii"                    
[365] "Rhacophorus_bengkuluensis"          "Rhacophorus_indonesiensis"         
[367] "Rhacophorus_malkmusi"               "Rhacophorus_rhyssocephalus"        
[369] "Rhinella_beebei"                    "Rhombophryne_botabota"             
[371] "Rhombophryne_longicrus"             "Rhombophryne_ornata"               
[373] "Rhombophryne_savaka"                "Rhombophryne_tany"                 
[375] "Rhombophryne_vaventy"               "Rohanixalus_baladika"              
[377] "Rohanixalus_hansenae"               "Rohanixalus_nauli"                 
[379] "Rohanixalus_vittatus"               "Sanguirana_acai"                   
[381] "Sanguirana_mearnsi"                 "Scaphiophryne_matsoko"             
[383] "Scinax_x.signatus"                  "Sclerophrys_arabica"               
[385] "Sclerophrys_buchneri"               "Sclerophrys_camerunensis"          
[387] "Sclerophrys_funerea"                "Sclerophrys_gracilipes"            
[389] "Sclerophrys_gutturalis"             "Sclerophrys_kassasii"              
[391] "Sclerophrys_latifrons"              "Sclerophrys_maculata"              
[393] "Sclerophrys_pentoni"                "Sclerophrys_pusilla"               
[395] "Sclerophrys_regularis"              "Sclerophrys_superciliaris"         
[397] "Sclerophrys_tihamica"               "Sclerophrys_togoensis"             
[399] "Sclerophrys_tuberosa"               "Sclerophrys_xeros"                 
[401] "Sigalegalephrynus_mandailinguensis" "Sigalegalephrynus_minangkabauensis"
[403] "Sphenophryne_allisoni"              "Sphenophryne_brevicrus"            
[405] "Sphenophryne_coggeri"               "Sphenophryne_crassa"               
[407] "Sphenophryne_dentata"               "Sphenophryne_magnitympanum"        
[409] "Sphenophryne_miniafia"              "Sphenophryne_rhododactyla"         
[411] "Sphenophryne_rubra"                 "Sphenophryne_schlaginhaufeni"      
[413] "Sphenophryne_similis"               "Sphenophryne_stenodactyla"         
[415] "Sphenophryne_thomsoni"              "Stumpffia_kibomena"                
[417] "Sylvirana_cubitalis"                "Sylvirana_guentheri"               
[419] "Sylvirana_maosonensis"              "Sylvirana_mortenseni"              
[421] "Sylvirana_nigrovittata"             "Sylvirana_roberti"                 
[423] "Sylvirana_spinulosa"                "Theloderma_vietnamense"            
[425] "Trachycephalus_typhonius"           "Uperodon_nagaoi"                   
[427] "Uperodon_obscurus"                  "Uperodon_palmatus"                 
[429] "Uperodon_taprobanicus"              "Uperodon_variegatus"               
[431] "Wijayarana_masonii"                 "Wijayarana_melasma"                
[433] "Wijayarana_modiglianii"             "Wijayarana_sumatrana"              
[435] "Xenopus_allofraseri"                "Xenopus_calcaratus"                
[437] "Xenopus_fischbergi"                 "Xenopus_mellotropicalis"           
[439] "Xenorhina_tillacki"                 "Zhangixalus_achantharrhena"        
[441] "Zhangixalus_amamiensis"             "Zhangixalus_arboreus"              
[443] "Zhangixalus_arvalis"                "Zhangixalus_aurantiventris"        
[445] "Zhangixalus_dennysi"                "Zhangixalus_dulitensis"            
[447] "Zhangixalus_moltrechti"             "Zhangixalus_owstoni"               
[449] "Zhangixalus_prasinatus"             "Zhangixalus_prominanus"            
[451] "Zhangixalus_schlegelii"             "Zhangixalus_taipeianus"            
[453] "Zhangixalus_viridis"                "Zhangixalus_yinggelingensis" 
