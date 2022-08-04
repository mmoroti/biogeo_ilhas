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
options(max.print=350)
match.phylo.comm(amphibia_phy,anura_islands)

# nomenclature adjustments 
# genus
{names(anura_islands) <- gsub("Hylomantis","Agalychnis", fixed=T, names(anura_islands)) 
names(anura_islands) <- gsub("Boana", "Hypsiboas", fixed=T, names(anura_islands)) 
names(anura_islands) <- gsub("Amnirana", "Hylarana", fixed=T, names(anura_islands)) 
names(anura_islands) <- gsub("Cophyla", "Platypelis", fixed=T, names(anura_islands)) 
names(anura_islands) <- gsub("Cornufer", "Platymantis", fixed=T, names(anura_islands)) 
names(anura_islands) <- gsub("Dryophytes", "Hyla", fixed=T, names(anura_islands)) 
names(anura_islands) <- gsub("Zhangixalus", "Rhacophorus", fixed=T, names(anura_islands)) 
names(anura_islands) <- gsub("Lithobates", "Rana", fixed=T, names(anura_islands)) 
names(anura_islands) <- gsub("Minervarya", "Fejervarya", fixed=T, names(anura_islands)) 
names(anura_islands) <- gsub("Chalcorana", "Hylarana", fixed=T, names(anura_islands)) 
names(anura_islands) <- gsub("Leptobrachella", "Leptolalax", fixed=T, names(anura_islands))
names(anura_islands) <- gsub("Ololygon", "Scinax", fixed=T, names(anura_islands))}
# epitete
{names(anura_islands) <- gsub("Hypsiboas_albomarginata", "Hypsiboas_albomarginatus", fixed=T, names(anura_islands)) 
names(anura_islands) <- gsub("Hypsiboas_calcarata", "Hypsiboas_calcaratus", fixed=T, names(anura_islands)) 
names(anura_islands) <- gsub("Hypsiboas_fasciata", "Hypsiboas_fasciatus", fixed=T, names(anura_islands)) 
names(anura_islands) <- gsub("Hypsiboas_geographica", "Hypsiboas_geographicus", fixed=T, names(anura_islands)) 
names(anura_islands) <- gsub("Hypsiboas_pulchella", "Hypsiboas_pulchellus", fixed=T, names(anura_islands)) 
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
names(anura_islands) <- gsub("Callimedusa_tomopterna","Phyllomedusa_tomopterna", fixed=T, names(anura_islands)) 
names(anura_islands) <- gsub("Chiasmocleis_hudsoni","Syncope_hudsoni", fixed=T, names(anura_islands))
names(anura_islands) <- gsub("Litoria_alboguttata", "Cyclorana_alboguttata", fixed=T, names(anura_islands))
names(anura_islands) <- gsub("Litoria_australis", "Cyclorana_australis", fixed=T, names(anura_islands))
names(anura_islands) <- gsub("Litoria_brevipes", "Cyclorana_brevipes", fixed=T, names(anura_islands))
names(anura_islands) <- gsub("Litoria_brevipes", "Cyclorana_brevipes", fixed=T, names(anura_islands))
names(anura_islands) <- gsub("Litoria_cryptotis", "Cyclorana_cryptotis", fixed=T, names(anura_islands))
names(anura_islands) <- gsub("Litoria_longipes", "Cyclorana_longipes", fixed=T, names(anura_islands))
names(anura_islands) <- gsub("Litoria_novaehollandiae", "Cyclorana_novaehollandiae", fixed=T, names(anura_islands))
names(anura_islands) <- gsub("Litoria_nullicedens", "Cyclorana_nullicedens", fixed=T, names(anura_islands))
names(anura_islands) <- gsub("Litoria_vagitus", "Cyclorana_vagitus", fixed=T, names(anura_islands))}

#Prune with phylogeny Pyron and Wiens consensus
phy_anura <- prune.sample(anura_islands, amphibia_phy) # 2159 spp
plot(phy_anura, type="fan",show.tip.label = FALSE)

#quantas espécies perdemos?
2503-2159 #344
#save phy islands
write.tree(phy_anura, file = "anura_islands.tre")

#"Dropping taxa from the community because they are not present in the phylogeny                
#[3] "Aglyptodactylus_australis"     "Aglyptodactylus_chorus"       
#[5] "Aglyptodactylus_inguinalis"    "Amolops_albispinus"           
#[7] "Anodonthyla_boulengeri"        "Ansonia_phuketensis"          
#[9] "Ansonia_teneritas"             "Ansonia_vidua"                
#[11] "Aphantophryne_nana"            "Asterophrys_eurydactyla"      
#[13] "Asterophrys_foja"              "Asterophrys_marani"           
#[15] "Asterophrys_pullifer"          "Asterophrys_slateri"          
#[17] "Austrochaperina_alexanderi"    "Austrochaperina_laurae"       
#[19] "Austrochaperina_punctata"      "Austrochaperina_rudolfarndti" 
#[21] "Bijurana_nicobariensis"        "Blommersia_nataliae"          
#[23] "Blommersia_transmarina"        "Boophis_ankarafensis"         
#[25] "Boophis_boppa"                 "Boophis_nauticus"             
#[27] "Buergeria_choui"               "Buergeria_otai"               
#[29] "Callulops_argus"               "Callulops_biakensis"          
#[31] "Callulops_bicolor"             "Callulops_neuhaussi"          
#[33] "Callulops_stellatus"           "Callulops_taxispilotus"       
#[35] "Callulops_wondiwoiensis"       "Callulops_yapenensis"         
#[37] "Hylarana_parvaccola"           "Chirixalus_nongkhorensis"     
#[39] "Chirixalus_trilaksonoi"        "Choerophryne_alainduboisi"    
#[41] "Choerophryne_alpestris"        "Choerophryne_bickfordi"       
#[43] "Choerophryne_bisyllaba"        "Choerophryne_brevicrus"       
#[45] "Choerophryne_brunhildae"       "Choerophryne_crucifer"        
#[47] "Choerophryne_darlingtoni"      "Choerophryne_epirrhina"       
#[49] "Choerophryne_exclamitans"      "Choerophryne_fafniri"         
#[51] "Choerophryne_grylloides"       "Choerophryne_gudrunae"        
#[53] "Choerophryne_gunnari"          "Choerophryne_laurini"         
#[55] "Choerophryne_multisyllaba"     "Choerophryne_murrita"         
#[57] "Choerophryne_pandanicola"      "Choerophryne_pipiens"         
#[59] "Choerophryne_rhenaurum"        "Choerophryne_sanguinopicta"   
#[61] "Choerophryne_siegfriedi"       "Choerophryne_swanhildae"      
#[63] "Choerophryne_tubercula"        "Choerophryne_valkuriarum"     
#[65] "Choerophryne_variegata"        "Cophixalus_cateae"            
#[67] "Cophixalus_hannahae"           "Cophixalus_rajampatensis"     
#[69] "Cophixalus_salawatiensis"      "Cophixalus_viridis"           
#[71] "Platypelis_berara"             "Platypelis_karenae"           
#[73] "Platypelis_maharipeo"          "Platypelis_noromalalae"       
#[75] "Platypelis_occultans"          "Platypelis_phyllodactyla"     
#[77] "Platypelis_puellarum"          "Platypelis_rava"              
#[79] "Copiula_alpestris"             "Copiula_annanoreenae"         
#[81] "Copiula_guttata"               "Copiula_lennarti"             
#[83] "Platymantis_bufoniformis"      "Platymantis_custos"           
#[85] "Platymantis_elegans"           "Platymantis_exedrus"          
#[87] "Platymantis_gigas"             "Platymantis_hedigeri"         
#[89] "Platymantis_heffernani"        "Platymantis_malukuna"         
#[91] "Platymantis_mediodiscus"       "Platymantis_minutus"          
#[93] "Platymantis_opisthodon"        "Platymantis_trossulus"        
#[95] "Platymantis_vertebralis"       "Platymantis_vogti"            
#[97] "Platymantis_wolfi"             "Dendropsophus_baileyi"        
#[99] "Hyla_cinereus"                 "Hyla_gratiosus"               
#[101] "Hyla_immaculatus"              "Hyla_japonicus"               
#[103] "Hyla_squirellus"               "Eleutherodactylus_campi"      
#[105] "Fejervarya_pulla"              "Glandirana_emeljanovi"        
#[107] "Glandirana_rugosa"             "Glandirana_susurra"           
#[109] "Glyphoglossus_brooksii"        "Glyphoglossus_capsus"         
#[111] "Glyphoglossus_flavus"          "Glyphoglossus_guttulatus"     
#[113] "Glyphoglossus_smithi"          "Glyphoglossus_volzi"          
#[115] "Guibemantis_diphonus"          "Hydrophylax_gracilis"         
#[117] "Hydrophylax_leptoglossa"       "Hydrophylax_malabaricus"      
#[119] "Hyla_perrini"                  "Hylophorbus_infulatus"        
#[121] "Hylophorbus_sigridae"          "Hyperolius_koehleri"          
#[123] "Hyperolius_stictus"            "Hypopachus_ustus"             
#[125] "Indosylvirana_aurantiaca"      "Indosylvirana_serendipi"      
#[127] "Indosylvirana_temporalis"      "Kalophrynus_barioensis"       
#[129] "Lankanectes_pera"              "Leptolalax_baluensis"         
#[131] "Leptolalax_brevicrus"          "Leptolalax_fuliginosa"        
#[133] "Leptolalax_itiokai"            "Leptolalax_juliandringi"      
#[135] "Leptolalax_laui"               "Leptolalax_marmorata"         
#[137] "Leptolalax_maura"              "Leptolalax_mjobergi"          
#[139] "Leptolalax_natunae"            "Leptolalax_palmata"           
#[141] "Leptolalax_parva"              "Leptolalax_picta"             
#[143] "Leptolalax_sabahmontana"       "Leptolalax_serasanae"         
#[145] "Leptobrachium_kantonishikawai" "Leptodactylus_pallidirostris" 
#[147] "Leptophryne_javanica"          "Limnonectes_cintalubang"      
#[149] "Limnonectes_larvaepartus"      "Limnonectes_quangninhensis"   
#[151] "Limnonectes_rhacodus"          "Rana_catesbeianus"            
#[153] "Rana_sphenocephalus"           "Rana_sylvaticus"              
#[155] "Litoria_alboguttata"           "Litoria_australis"            
#[157] "Litoria_brevipes"              "Litoria_cryptotis"            
#[159] "Litoria_longipes"              "Litoria_novaehollandiae"      
#[161] "Litoria_nullicedens"           "Litoria_pallidofemora"        
#[163] "Litoria_pinocchio"             "Litoria_rueppelli"            
#[165] "Litoria_vagitus"               "Liuixalus_calcarius"          
#[167] "Lysapsus_bolivianus"           "Mantophryne_insignis"         
#[169] "Mantophryne_menziesi"          "Megophrys_aceras"             
#[171] "Megophrys_baluensis"           "Megophrys_boettgeri"          
#[173] "Megophrys_brachykolos"         "Megophrys_carinense"          
#[175] "Megophrys_dringi"              "Megophrys_edwardinae"         
#[177] "Megophrys_major"               "Megophrys_microstoma"         
#[179] "Megophrys_parallela"           "Megophrys_parva"              
#[181] "Meristogenys_penrissenensis"   "Microhyla_kuramotoi"          
#[183] "Microhyla_mihintalei"          "Micryletta_erythropoda"       
#[185] "Micryletta_nigromaculata"      "Fejervarya_chiangmaiensis"    
#[187] "Fejervarya_muangkanensis"      "Nanohyla_perparva"            
#[189] "Nanohyla_petrigena"            "Nidirana_adenopleura"         
#[191] "Nidirana_hainanensis"          "Nidirana_okinavana"           
#[193] "Nyctimystes_infrafrenatus"     "Nyctimystes_myolae"           
#[195] "Nyctimystes_ocreptus"          "Nyctimystes_sauroni"          
#[197] "Nyctimystes_tyleri"            "Ololygon_agilis"              
#[199] "Ololygon_argyreornata"         "Ololygon_berthae"             
#[201] "Ololygon_brieni"               "Ololygon_catharinae"          
#[203] "Ololygon_flavoguttata"         "Ololygon_humilis"             
#[205] "Ololygon_jureia"               "Ololygon_littoralis"          
#[207] "Ololygon_perpusilla"           "Ololygon_rizibilis"           
#[209] "Ololygon_trapicheiroi"         "Oreophryne_albitympanum"      
#[211] "Oreophryne_albomaculata"       "Oreophryne_anser"             
#[213] "Oreophryne_aurora"             "Oreophryne_banshee"           
#[215] "Oreophryne_brunnea"            "Oreophryne_cameroni"          
#[217] "Oreophryne_choerophrynoides"   "Oreophryne_curator"           
#[219] "Oreophryne_equus"              "Oreophryne_flavomaculata"     
#[221] "Oreophryne_gagneorum"          "Oreophryne_graminis"          
#[223] "Oreophryne_lemur"              "Oreophryne_matawan"           
#[225] "Oreophryne_meliades"           "Oreophryne_nicolasi"          
#[227] "Oreophryne_parkopanorum"       "Oreophryne_penelopeia"        
#[229] "Oreophryne_philosylleptoris"   "Oreophryne_phoebe"            
#[231] "Oreophryne_picticrus"          "Oreophryne_pseudunicolor"     
#[233] "Oreophryne_roedeli"            "Paedophryne_titan"            
#[235] "Papurana_arfaki"               "Papurana_attigua"             
#[237] "Papurana_aurata"               "Papurana_daemeli"             
#[239] "Papurana_elberti"              "Papurana_florensis"           
#[241] "Papurana_garritor"             "Papurana_grisea"              
#[243] "Papurana_jimiensis"            "Papurana_kreffti"             
#[245] "Papurana_milleti"              "Papurana_milneana"            
#[247] "Papurana_moluccana"            "Papurana_novaeguineae"        
#[249] "Papurana_papua"                "Papurana_supragrisea"
#[251] "Papurana_volkerjane"                "Papurana_waliesa"             
#[253] "Pelobates_balcanicus"               "Pelophryne_penrissenensis"         
#[255] "Petropedetes_newtonii"              "Philautus_catbaensis"              
#[257] "Philautus_nephophilus"              "Phrynoidis_asper"                  
#[259] "Pithecopus_hypochondrialis"         "Platymantis_quezoni"               
#[261] "Plethodontohyla_alluaudi"           "Polypedates_pseudotilophus"        
#[263] "Pseudophilautus_conniffae"          "Pulchrana_baramica"                
#[265] "Pulchrana_centropeninsularis"       "Pulchrana_debussyi"                
#[267] "Pulchrana_glandulosa"               "Pulchrana_grandocula"              
#[269] "Pulchrana_guttmani"                 "Pulchrana_laterimaculata"          
#[271] "Pulchrana_mangyanum"                "Pulchrana_melanomenta"             
#[273] "Pulchrana_moellendorffi"            "Pulchrana_picturata"               
#[275] "Pulchrana_rawa"                     "Pulchrana_siberu"                  
#[277] "Pulchrana_signata"                  "Pulchrana_similis"                 
#[279] "Rana_neba"                          "Rentapia_everetti"                 
#[281] "Rentapia_hosii"                     "Rhacophorus_bengkuluensis"         
#[283] "Rhacophorus_indonesiensis"          "Rhacophorus_malkmusi"              
#[285] "Rhacophorus_rhyssocephalus"         "Rhinella_beebei"                   
#[287] "Rhombophryne_botabota"              "Rhombophryne_longicrus"            
#[289] "Rhombophryne_ornata"                "Rhombophryne_savaka"               
#[291] "Rhombophryne_tany"                  "Rhombophryne_vaventy"              
#[293] "Rohanixalus_baladika"               "Rohanixalus_hansenae"              
#[295] "Rohanixalus_nauli"                  "Rohanixalus_vittatus"              
#[297] "Sanguirana_acai"                    "Sanguirana_mearnsi"                
#[299] "Scaphiophryne_matsoko"              "Scinax_x.signatus"                 
#[301] "Sclerophrys_arabica"                "Sclerophrys_buchneri"              
#[303] "Sclerophrys_camerunensis"           "Sclerophrys_funerea"               
#[305] "Sclerophrys_gracilipes"             "Sclerophrys_gutturalis"            
#[307] "Sclerophrys_kassasii"               "Sclerophrys_latifrons"             
#[309] "Sclerophrys_maculata"               "Sclerophrys_pentoni"               
#[311] "Sclerophrys_pusilla"                "Sclerophrys_regularis"             
#[313] "Sclerophrys_superciliaris"          "Sclerophrys_tihamica"              
#[315] "Sclerophrys_togoensis"              "Sclerophrys_tuberosa"            
#[317] "Sclerophrys_xeros"                "Sigalegalephrynus_mandailinguensis"
#[319] "Sigalegalephrynus_minangkabauensis" "Sphenophryne_allisoni"             
#[321] "Sphenophryne_brevicrus"             "Sphenophryne_coggeri"              
#[323] "Sphenophryne_crassa"                "Sphenophryne_dentata"              
#[325] "Sphenophryne_magnitympanum"         "Sphenophryne_miniafia"             
#[327] "Sphenophryne_rhododactyla"          "Sphenophryne_rubra"                
#[329] "Sphenophryne_schlaginhaufeni"       "Sphenophryne_similis"              
#[331] "Sphenophryne_stenodactyla"          "Sphenophryne_thomsoni"             
#[333] "Stumpffia_kibomena"                 "Sylvirana_cubitalis"               
#[335] "Sylvirana_guentheri"                "Sylvirana_maosonensis"             
#[337] "Sylvirana_mortenseni"               "Sylvirana_nigrovittata"            
#[339] "Sylvirana_roberti"                  "Sylvirana_spinulosa"               
#[341] "Theloderma_vietnamense"             "Trachycephalus_typhonius"          
#[343] "Uperodon_nagaoi"                    "Uperodon_obscurus"                 
#[345] "Uperodon_palmatus"                  "Uperodon_taprobanicus"             
#[347] "Uperodon_variegatus"                "Wijayarana_masonii"                
#[349] "Wijayarana_melasma"                 "Wijayarana_modiglianii"            
#[351] "Wijayarana_sumatrana"               "Xenopus_allofraseri"               
#[353] "Xenopus_calcaratus"                 "Xenopus_fischbergi"                
#[355] "Xenopus_mellotropicalis"            "Xenorhina_tillacki"                
#[357] "Rhacophorus_amamiensis"  
