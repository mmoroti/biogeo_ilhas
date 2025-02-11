# Tue Feb 11 15:36:08 2025 ------------------------------
# Packages ----
library(FD)
library(tidyverse)
library(picante)
library(tidyverse)
library(rgbif)
library(phyloregion)
library(gawdis)
library(funspace)
library(maditr)
library(ggExtra)

# Set local directory before
local_directory <- file.path(
  "E:", "datasets_centrais")

# Datasets ----
# TetrapodTraits
tetrapod_traits <-  read.csv2(
  file.path(local_directory, "TetrapodTraits_1.0.0.csv"), sep = ",") 

amphibio_traits <- read.csv2(
  file.path(local_directory, "AmphiBIO_v1",
            "AmphiBIO_v1.csv"), sep = ",")

# Phylogeny of amphibians (Jetz and Pyron, 2018)
amphibia_tree <- read.tree(file.path(
  local_directory, "amph_shl_new_Consensus_7238.tre"))

# Community of IUCN polygons
polygons_data <- as_tibble(read.csv2(
  file.path(local_directory, "biogeo_islands",
            "anura_islands_without0.csv"), header = TRUE, sep = ";"))[,-1]

# Community point of occurences data from GBIF
occurences_data <- read_delim(
  file.path(local_directory, "biogeo_islands",
            "data_gbif.csv"), delim = ",")[,-1]

# Species list Tropical and Temperate domains
amphibians_world <- as_tibble(
  read.csv2(file.path(
    local_directory, "biogeo_islands",
    "cross_validation", "anura_func_new.csv"
  ), sep=",")) %>%
  select(
    "binomial", "climate"
  )

# Point of occurences in islands 
list_species_gbif <- read_delim(
  file.path(local_directory, "biogeo_islands",
            "cross_validation", "data_gbif.csv"),
  delim = ",")[,-1]

# Nomenclature harmonize ----
# First, we obtain the speciesKey for each specie in each of our dataset
# Phylogeny species names
names_phy <- amphibia_tree$tip.label %>%
  name_backbone_checklist()

names_phy_adj <- names_phy %>%
  select(verbatim_name, speciesKey) %>%
  rename(names_phy = verbatim_name)

# Community species 
# polygon tem 1924 spp
names_community <- names(polygons_data) %>%
  name_backbone_checklist()

# Traits AmphiBIO
names_amphibio <- amphibio_traits$Species %>%
  name_backbone_checklist()

names_amphibio_adj <- names_amphibio %>%
  select(verbatim_name, speciesKey) %>%
  rename(Species = verbatim_name) %>%
  left_join(amphibio_traits, by = "Species")

# Amphibia tropical vs climate 
species_list_islands <- amphibians_world$binomial %>%
  name_backbone_checklist()

names_islands_adj <- species_list_islands %>%
  select(verbatim_name, speciesKey) %>%
  rename(Species = verbatim_name)

# GBIF id
list_species_gbif_key <- list_species_gbif %>%
  distinct(species, .keep_all = TRUE) %>%
  pull(species) %>%
  name_backbone_checklist()

list_gbif <- list_species_gbif %>%
  distinct(species, .keep_all = TRUE) %>%
  select(species, family) %>%
  # filter gymnophiona and caudata
  filter(family != "Ambystomatidae",
         family != "Hynobiidae", 
         family != "Salamandridae", 
         family != "Dermophiidae", 
         family != "Plethodontidae", 
         family != "Indotyphlidae",
         family != "Ichthyophiidae",
         family != "Cryptobranchidae",
         family != "Proteidae",
         family != "Amphiumidae",
         family != "Siphonopidae"
  ) 

list_gbif_key <- list_gbif %>%  
  left_join(list_species_gbif_key, by = c("species" = "verbatim_name")) %>%
  select(species, speciesKey) 

### Communities (Polygons) ----
# vamos selecionar apenas os specieskey nao repetidos
# According Frost, 2025, these 3 duplicates speciesKey are synonyms of each other
# Chalcorana kampeni is Chalcorana_crassiovis
# Nyctimystes_oktediensis is Nyctimystes_disruptus
# Peltophryne_fracta is Peltophryne_guentheri
synonyms <- c("Chalcorana_kampeni", "Nyctimystes_oktediensis",  "Peltophryne_fracta")

names_community_adj <- names_community %>%
  select(scientificName, verbatim_name, speciesKey) %>%
  rename(names_community = verbatim_name) %>%
  filter(!names_community %in% synonyms) %>%
  filter(names_community != "id") # 1921 species

# sem duplicatas nas speciesKey
names_community_adj %>%
  filter(duplicated(speciesKey) | duplicated(speciesKey, fromLast = TRUE)) 

# harmonizando filogenia & comunidades
# esses nomes sao os nomes respectivos que queremos na filogenia, no
# trabalho do pyron tem varias especies consideradas sinonimos pelo 
# backbone do gbif, isso eh para garantir que estamos acessando o mesmo nome
# que temos na comunidade e na
names_phy_select <- c(
  "Adenomus_kandianus",
  "Hylarana_albolabris",
  "Arthroleptis_stenodactylus",
  "Aubria_occidentalis",
  "Bufo_gargarizans",
  "Hylarana_crassiovis",
  "Platymantis_boulengeri",
  "Hyla_japonica",
  "Espadarana_prosoblepon",
  "Eupsophus_calcaratus",
  "Fejervarya_cancrivora",
  "Hyperolius_glandicolor",
  "Kaloula_pulchra",
  "Kurixalus_eiffingeri",
  "Leptopelis_argenteus",
  "Leptopelis_flavomaculatus",
  "Limnonectes_laticeps",
  "Microhyla_berdmorei",
  "Babina_adenopleura",
  "Nyctimystes_disruptus",
  "Pelophryne_guentheri",
  "Phrynobatrachus_latifrons",
  "Platymantis_guentheri",
  "Ptychadena_nilotica",
  "Pedostibes_everetti",
  "Rhacophorus_annamensis",
  "Rhinella_crucifer",
  "Sanguirana_sanguinea",
  "Amietophrynus_togoensis",
  "Hylarana_nigrovittata",
  "Teratohyla_pulverata",
  "Uperoleia_rugosa",
  "Peltophryne_guentheri"
) 

# retirar duplicatas
rename_cols <- left_join(names_community_adj,
                         names_phy_adj,
                         by = "speciesKey") %>%
  group_by(speciesKey) %>%
  filter(n() == 1) %>% # Mantém apenas speciesKey que aparecem uma vez
  ungroup()

# filtrar para as especies selecionadas nas arvores e depois readicionar
# ao conjunto de dados de cima!
duplicates <- left_join(names_community_adj,
                        names_phy_adj,
                        by = "speciesKey") %>%
  filter(duplicated(speciesKey) | duplicated(speciesKey, fromLast = TRUE)) %>%
  filter(names_phy %in% names_phy_select) 

community_names_with_speciesKey <- bind_rows(
  rename_cols,
  duplicates
)

community_adj_names <- bind_rows(
  rename_cols,
  duplicates
) %>%
  select(names_phy, names_community) %>%
  remove_missing() %>% # sp nao encontradas nas arvores
  deframe()

# ok
# 1918 nomes identificados na arvore, representacao de quase 99%
# Sem posicao na filogenia
remove_species <- c("Bijurana_nicobariensis",
                    "Dryophytes_femoralis", 
                    "Lithobates_palustris")

# Unifying synonyms 
polygons_data$Chalcorana_crassiovis <- pmax(community_data$Chalcorana_kampeni,
                                            community_data$Chalcorana_crassiovis)

polygons_data$Nyctimystes_disruptus <- pmax(community_data$Nyctimystes_oktediensis,
                                            community_data$Nyctimystes_disruptus)

polygons_data$Peltophryne_guentheri <- pmax(community_data$Peltophryne_fracta,
                                            community_data$Peltophryne_guentheri)

polygons_adj <- polygons_data %>%
  select(-all_of(c(synonyms, remove_species))) %>%
  rename(all_of(community_adj_names)) 

# ncol(polygons_adj) # 1918 

### Communities (Points of occurence) ---- 
# Unifying GBIF and Traits
list_gbif_phy <- left_join(list_gbif_key, 
                           names_phy_adj,
                           by = 'speciesKey') %>%
  remove_missing() # 76 spp sem posicao na filogenia e sem traits

synonyms_gbif <- c("Hyla_japonica",
                   "Bufo_gargarizans",
                   "Platymantis_guentheri",
                   "Rhinella_crucifer",
                   "Fejervarya_cancrivora",
                   "Kaloula_pulchra",
                   "Nyctimystes_disruptus",
                   "Pelophylax_lessonae",
                   "Limnonectes_khasianus",
                   "Babina_adenopleura",
                   "Hylarana_crassiovis",
                   "Sanguirana_sanguinea",
                   "Ptychadena_nilotica",
                   "Kurixalus_eiffingeri",
                   "Hyla_immaculata",
                   "Hylarana_albolabris",
                   "Arthroleptis_stenodactylus",
                   "Platymantis_boulengeri",
                   "Aubria_subsigillata",
                   "Aubria_subsigillata",
                   "Microhyla_berdmorei",
                   "Neobatrachus_sudelli",
                   "Odorrana_tiannanensis",
                   "Sphaenorhynchus_platycephalus",
                   "Hylarana_nigrovittata",
                   "Pedostibes_everetti")

list_gbif_phy %>%
  filter(duplicated(speciesKey) | duplicated(speciesKey, fromLast = TRUE)) %>%
  View()

# retirar duplicatas
rename_cols <- list_gbif_phy %>%
  group_by(speciesKey) %>%
  filter(n() == 1) %>% # Mantém apenas speciesKey que aparecem uma vez
  ungroup()

# filtrar para as especies selecionadas nas arvores e depois readicionar
# ao conjunto de dados de cima!
duplicates <- list_gbif_phy %>%
  filter(duplicated(speciesKey) | duplicated(speciesKey, fromLast = TRUE)) %>%
  filter(names_phy %in% synonyms_gbif) 

community_gbif_with_speciesKey <- bind_rows(
  rename_cols,
  duplicates
) 

# transform em communities
matrix_species_gbif <- as.data.frame(dcast(
  list_species_gbif, OBJECTID_1 ~ species, length)) %>%
  filter(OBJECTID_1 != "17831",
         OBJECTID_1 != "2920",
         OBJECTID_1 != "2940",
         OBJECTID_1 != "2959",
         OBJECTID_1 != "3278",
         OBJECTID_1 != "3372",
         OBJECTID_1 != "17831",
         OBJECTID_1 != "5017") %>% # remove island with 0 spp.
  column_to_rownames("OBJECTID_1") %>% 
  mutate(across(everything(), ~ ifelse(. > 0, 1, .))) 

matrix_gbif <- matrix_species_gbif %>% 
  select(community_gbif_with_speciesKey$species)

# retirando ilhas com soma 0
matrix_gbif_without0 <- matrix_gbif[rowSums(matrix_gbif) >= 1, ]

# vetor para renomear seguindo a filogenia
list_rename_gbif <- community_gbif_with_speciesKey %>%
  select(names_phy, species) %>%
  deframe()
# renomeando segundo a filogenia
matrix_gbif_renamed <- matrix_gbif_without0 %>%
  rename(all_of(list_rename_gbif)) %>%
  select(order(names(.))) 

### AmphiBIO ----
traits_amphibio <- names_amphibio_adj %>%
  select(Species, speciesKey, Dir, Lar, Viv)

# esses nomes sao os nomes respectivos que queremos no AmphiBIO, no
# AmphiBIO também tem varias especies consideradas sinonimos pelo 
# backbone do gbif, isso eh para garantir que estamos acessando o mesmo nome
# que temos na comunidade e no AmphiBIO
names_amp_select <- c(
  "Cophixalus ateles",
  "Cophixalus verrucosus",
  "Microhyla borneensis",
  "Rana sauteri",
  "Adenomus kandianus",
  "Hylarana albolabris",
  "Bufo gargarizans",
  "Hylarana crassiovis",
  "Hyla japonica",
  "Leptopelis flavomaculatus",
  "Limnonectes laticeps",
  "Microhyla berdmorei",
  "Babina adenopleura",
  "Litoria disrupta",
  "Peltophryne guentheri",
  "Phrynobatrachus latifrons",
  "Pedostibes everetti",
  "Rhinella crucifer",
  "Sanguirana sanguinea",
  "Amietophrynus togoensis",
  "Hylarana nigrovittata",
  "Uperoleia rugosa"
) 

# retirar duplicatas
species_amphibio <- left_join(community_names_with_speciesKey,
                              traits_amphibio,
                              by = "speciesKey") %>%
  group_by(speciesKey) %>%
  filter(n() == 1) %>% # Mantém apenas speciesKey que aparecem uma vez
  ungroup()

# filtrar para as especies selecionadas nas arvores e depois readicionar
# ao conjunto de dados de cima!
duplicates_amphibio <- left_join(community_names_with_speciesKey,
                                 traits_amphibio,
                                 by = "speciesKey") %>%
  filter(duplicated(speciesKey) | duplicated(speciesKey, fromLast = TRUE)) %>%
  filter(Species %in% names_amp_select)  

amphibio_names_with_speciesKey <- bind_rows(
  species_amphibio,
  duplicates_amphibio
) %>%
  filter(!names_phy %in% remove_species) %>% # species without phylogeny
  select(names_phy, Dir, Lar, Viv) %>%
  mutate(Scientific.Name = gsub("_", " ", names_phy))

### TetrapodTraits ----
list_polygons <- polygons_adj %>%
  gather(key="Scientific.Name",value="count",
         Hylarana_luctuosa:Rhacophorus_viridis) %>% 
  distinct(Scientific.Name, .keep_all = TRUE) %>%
  select(-id, -count) %>%
  mutate(Scientific.Name = gsub("_", " ", Scientific.Name))

names(tetrapod_traits)

tetrapod_traits_select <- tetrapod_traits %>%
  filter(Class == "Amphibia") %>%
  select("Scientific.Name", "IUCN_Binomial",
         "BodyLength_mm", "ImputedLength",
         "Diu", "Noc", "ImputedActTime", 
         "Fos", "Ter", "Aqu", "Arb", "ImputedHabitat")


### Species list tropical, temperate & both ----
without_duplicates <- amphibians_world %>%
  distinct() 

names_islands_dup <- names_islands_adj %>%
  distinct() 

View(names_islands_dup)

# identificando as especies pelo speciesKey
species_climate_islands <- left_join(
  names_islands_dup,
  without_duplicates, by = c("Species" = "binomial"))

# identificando quais especies estao em cada lugar
communities_islands <- left_join(community_adj_names,
                                 species_climate_islands, by = "speciesKey") %>%
  select(names_phy, climate) 

# adicionando informacao da posicao da filha para as especies com traits
# 1918 spp. Mas algumas ocorrem nos dois lugares, por isso o left_join ira
# duplicar algumas especies
trait_list <- data.frame(Species = row.names(amp_islands_traits))

species_per_islands <- left_join(trait_list,
                                 communities_islands, by = c("Species" = "names_phy"))

# algumas especies tem nos dois (temperada e tropical)
#species_per_islands %>%
#  distinct(Species, .keep_all = TRUE) %>% View()

species_per_islands_unique <-species_per_islands %>%
  filter(!duplicated(Species) & !duplicated(Species, fromLast = TRUE)) 

species_per_islands_both <- species_per_islands %>%
  filter(duplicated(Species) | duplicated(Species, fromLast = TRUE)) %>%
  mutate(climate = "both") %>%
  distinct(Species, .keep_all = TRUE)

species_per_islands_classified <- bind_rows(
  species_per_islands_unique,
  species_per_islands_both
)


# Metrics ----
### Obtain phylogenetic diversity ----
match.phylo.comm(amphibia_tree, polygons_adj)

polygons_phy <- prune.sample(polygons_adj, amphibia_tree) 

# Phylogenetic diversity
# filter dataframe
matrix <- polygons_adj %>%
  column_to_rownames("id") %>%
  select(sort(names(.)))

# alguma comunidade com zero?
any(rowSums(matrix) == 0)

# Analysis of phylogenetic endemism and phyloregions 
#Sparse community matrix - requirement from analysis of phyloregion package
matrix_sparse_iucn <- dense2sparse(matrix)

# phylo_endemism
pe_iucn <- phylo_endemism(matrix_sparse_iucn, polygons_phy)
head(data.frame(pe_iucn))

### Obtain functional diversity ----
amphibia_islands_traits <- left_join(list_polygons, 
                                     tetrapod_traits_select, by = 'Scientific.Name') %>%
  mutate(across(c(Diu, Noc, Fos, Ter, Aqu, Arb), ~ ifelse(. > 0.5, 1, 0))) 

amp_islands_traits_complete <- left_join(amphibia_islands_traits,
                                         amphibio_names_with_speciesKey,
                                         by = "Scientific.Name") %>%
  arrange(Scientific.Name)

# We need to complete species without development mode
# write.csv2(amp_islands_traits_complete, "amp_islands_traits_tofill.csv")

# Load dataset complete
amp_islands_traits_complete <- read.csv2(
  file.path(local_directory, "biogeo_islands",
            "amp_islands_traits_complete.csv")) %>%
  mutate(Scientific.Name = gsub(" ", "_", Scientific.Name)) %>%
  arrange(Scientific.Name)

columns_flags <- c("IUCN_Binomial", "Ref_Mode_Reproductive", "ImputedHabitat",
                   "ImputedActTime", "ImputedLength", "Infered")

# With imputation
amp_islands_traits <- amp_islands_traits_complete %>%
  select(-all_of(columns_flags)) %>%
  mutate(BodyLength_mm = log(BodyLength_mm)) %>%
  column_to_rownames("Scientific.Name") 
visdat::vis_miss(amp_islands_traits)

# Without imputation
amphibia_without_imputed <- amp_islands_traits_complete %>%
  filter(ImputedLength != 1 & ImputedActTime != 1 & 
           ImputedHabitat != 1 & Infered != 1) %>%
  select(-all_of(columns_flags)) %>%
  mutate(BodyLength_mm = log(BodyLength_mm)) %>%
  column_to_rownames("Scientific.Name") 

# gawdis distance
gaw.anura <-gawdis(amp_islands_traits,
                   w.type="analytic",
                   groups = c(1,
                              2,2,
                              3,3,3,3,
                              4,4,4))
# Community matrix
polygons_community <- polygons_adj %>%
  select(-id) %>%
  select(sort(names(.))) %>%
  as.matrix()

# FDis
disp_func_amphibia <- fdisp(
  gaw.anura, 
  polygons_community,
  tol = 1e-07)

functional_diversity <- dbFD(gaw.anura, polygons_community)

rich <- polygons_adj %>%
  column_to_rownames("id") %>%
  rowSums()

diversity_polygons <- data.frame(id= polygons_adj$id,
                                 FDisp = disp_func_amphibia$FDis,
                                 #FRich = functional_diversity$FRic,
                                 PhyEnd = pe_iucn,
                                 Richness = rich)
head(diversity_polygons)

### Obtain Trait-space with Funspace Package ----

# PCoA
trait_distance_gawdis_pcoa <- capscale(gaw.anura ~ 1,
                                       amp_islands_traits,
                                       distance = "gower")

# Autovalores (importância de cada eixo)
eigenvalues <- trait_distance_gawdis_pcoa$CA$eig

# Percentual de variância explicada por eixo
variance_explained <- eigenvalues / sum(eigenvalues)

# Coordenadas dos pontos nos eixos principais
ordination_scores <- scores(trait_distance_gawdis_pcoa, display = "sites")

# Calculando correlação entre os traços e os eixos da PCoA
trait_correlation <- cor(amp_islands_traits, ordination_scores, use = "pairwise.complete.obs")
print(trait_correlation)


# Funspace
trait_space_global <- funspace(x=trait_distance_gawdis_pcoa,
                               PCs=c(1,2),
                               n_divisions=1000)
fit <- envfit(ord = trait_distance_gawdis_pcoa,
              env = amp_islands_traits)

global_islands <- plot(trait_space_global,
                       type="global", # plot the global TPD
                       quant.plot=TRUE, # add quantile lines
                       arrows=TRUE, # add arrows for PCA loadings
                       arrows.length=0.9,
                       globalContour = T,
                       pnt = T,
                       pnt.col = rgb(0.8, 0.7, 0.1, alpha = 0.2),
                       axis.title.x = "Development mode (26.46%)",
                       axis.title.y = "Circadian Activity (22.01%)")
plot(fit, add = TRUE, col = 'black')

# adjust GAM
#y <- abs(scores(trait_distance_gawdis_pcoa, display = "sites")[, 1] * 
#           scores(trait_distance_gawdis_pcoa, display = "sites")[, 2] + 
#           rnorm(nrow(amp_islands_traits), 0, 1))

#fit.gam <- funspaceGAM(y = y, funspace = trait_space_global)

#plot(x = fit.gam, 
#     type = "global",
#     quant.plot = TRUE,
#     quant.col = "grey80")

#plot(fit, add = TRUE, col = 'black')

# islands position (temperate, tropical, both)
trait_space_climate <- funspace(
  x=trait_distance_gawdis_pcoa, 
  PCs=c(1,2),
  group.vec=species_per_islands_classified$climate,
  n_divisions=1000) 

# funspace per islands
plot(x=trait_space_climate, # funspace object
     type="groups", # plot the global TPD
     which.group = "temperate",
     quant.plot=TRUE, # add quantile lines
     arrows=TRUE, # add arrows for PCA loadings
     arrows.length=0.9,
     globalContour = T,
     pnt = T,
     pnt.col = rgb(0.8, 0.7, 0.1, alpha = 0.2),
     axis.title.x = "Development mode (26.46%)",
     axis.title.y = "Circadian Activity (22.01%)")
plot(fit, add = TRUE, col = 'black')

### 
tiff("global_funspace.tiff", width = 1700, height = 1700, res = 300)
plot(trait_space_global,
     type="global", # plot the global TPD
     quant.plot=TRUE, # add quantile lines
     arrows=TRUE, # add arrows for PCA loadings
     arrows.length=0.9,
     globalContour = T,
     pnt = T,
     pnt.col = rgb(0.8, 0.7, 0.1, alpha = 0.2),
     axis.title.line = 2.5,
     axis.title.x = "Development mode (26.46%)",
     axis.title.y = "Circadian Activity (22.01%)")
plot(fit, add = TRUE, col = 'black')
dev.off()

tiff("both_funspace.tiff", width = 1700, height = 1700, res = 300)
plot(x=trait_space_climate, # funspace object
     type="groups", # plot the global TPD
     which.group = "both",
     quant.plot = TRUE, # add quantile lines
     arrows = TRUE, # add arrows for PCA loadings
     arrows.length = 0.9,
     globalContour = TRUE,
     pnt = TRUE,
     pnt.col = rgb(0.8, 0.7, 0.1, alpha = 0.2),
     axis.title.line = 2.5,
     axis.title.x = "Development mode (26.46%)",
     axis.title.y = "Circadian Activity (22.01%)")
plot(fit, add = TRUE, col = 'black')
dev.off()

tiff("temperate_funspace.tiff", width = 1700, height = 1700, res = 300)
plot(x=trait_space_climate, # funspace object
     type="groups", # plot the global TPD
     which.group = "temperate",
     quant.plot = TRUE, # add quantile lines
     arrows = TRUE, # add arrows for PCA loadings
     arrows.length = 0.9,
     globalContour = TRUE,
     pnt = TRUE,
     pnt.col = rgb(0.8, 0.7, 0.1, alpha = 0.2),
     axis.title.line = 2.5,
     axis.title.x = "Development mode (26.46%)",
     axis.title.y = "Circadian Activity (22.01%)")
plot(fit, add = TRUE, col = 'black')
dev.off()

tiff("tropical_funspace.tiff", width = 1700, height = 1700, res = 300)
plot(x=trait_space_climate, # funspace object
     type="groups", # plot the global TPD
     which.group = "tropical",
     quant.plot = TRUE, # add quantile lines
     arrows = TRUE, # add arrows for PCA loadings
     arrows.length = 0.9,
     globalContour = TRUE,
     pnt = TRUE,
     pnt.col = rgb(0.8, 0.7, 0.1, alpha = 0.2),
     axis.title.line = 2.5,
     axis.title.x = "Development mode (26.46%)",
     axis.title.y = "Circadian Activity (22.01%)")
plot(fit, add = TRUE, col = 'black')
dev.off()

# funspace gam adjust
#fit.gam <- funspaceGAM(y = y, funspace = trait_space_climate)

#plot(x=fit.gam, # funspace object
#     type="groups", # plot the global TPD
#     quant.plot=TRUE, # add quantile lines
#     arrows=TRUE, # add arrows for PCA loadings
#     arrows.length=0.9) 

#plot(fit, add = TRUE, col = 'black')


# Cross-validation ----
remove <- colnames(matrix_gbif_renamed)[colSums(matrix_gbif_renamed) == 0]

amp_gbif_community <- matrix_gbif_renamed %>%
  select(-all_of(remove))

#### Evolutionary distincteness (Phylogenetic endemism) ----
gbif_phy <- prune.sample(amp_gbif_community, amphibia_tree) 
# alguma comunidade com zero?
any(rowSums(amp_gbif_community) == 0)
# Analysis of phylogenetic endemism and phyloregions 
#Sparse community matrix - requirement from analysis of phyloregion package
matrix_sparse_gbif <- dense2sparse(amp_gbif_community)
# phylo_endemism
pe_gbif <- phylo_endemism(matrix_sparse_gbif, gbif_phy)
# View(data.frame(pe_gbif))

### Functional diversity (Functional dispersion) ---- 
list_gbif_renamed <- data.frame('Scientific.Name' = names(amp_gbif_community))

tetrapod_traits_adj <- tetrapod_traits_select %>%
  mutate(Scientific.Name = gsub(" ", "_", Scientific.Name))  

amphibia_gbif_traits <- left_join(list_gbif_renamed, 
                                  tetrapod_traits_adj, by = 'Scientific.Name') %>%
  mutate(across(c(Diu, Noc, Fos, Ter, Aqu, Arb), ~ ifelse(. > 0.5, 1, 0)))

amp_list_traits <- amp_islands_traits_complete %>%
  select(Scientific.Name, Dir, Lar, Viv) 

amphibia_gbif_traits_tofill <- amphibia_gbif_traits %>%
  left_join(amp_list_traits, by = "Scientific.Name")

write.csv2(amphibia_gbif_traits_tofill, "amphibia_gbif_traits_tofill.csv" )

# Functional dispersion
amp_gbis_traits_complete <- read.csv2(
  file.path(local_directory, "biogeo_islands",
            "amphibia_gbif_traits_complete.csv")) %>%
  mutate(Scientific.Name = gsub(" ", "_", Scientific.Name),
         BodyLength_mm = log(as.numeric(BodyLength_mm))) %>%
  arrange(Scientific.Name) %>%
  select(-"IUCN_Binomial", -"X", -"Ref", -"ImputedHabitat",
         -"ImputedActTime", -"ImputedLength") 

nrow(amp_gbis_traits_complete) # 1494
ncol(matrix_gbif_renamed) # 1494 especies

amp_gbif_traits <- amp_gbis_traits_complete %>%
  filter(!Scientific.Name %in% remove) %>%
  arrange(Scientific.Name) %>%
  column_to_rownames("Scientific.Name")

nrow(amp_gbif_traits) # 1492
ncol(amp_gbif_community) # 1492 especies

# fazer a matriz de distancia
# gawdis distance
gaw.anura.gbif <-gawdis(amp_gbif_traits,
                        w.type="analytic",
                        groups = c(1,
                                   2,2,
                                   3,3,3,3,
                                   4,4,4))
# Community matrix
gbif_community <- amp_gbif_community %>%
  select(sort(names(.))) %>%
  as.matrix()

# calcular FDis
disp_func_amphibia_gbif <- fdisp(
  gaw.anura.gbif, 
  gbif_community,
  tol = 1e-07)

# save metrics validation
rich_gbif <- gbif_community %>%
  rowSums()

diversity_gbif <- data.frame(id= as.integer(rownames(gbif_community)),
                             FDisp_GBIF = disp_func_amphibia_gbif$FDis,
                             PhyEnd_GBIF = pe_gbif,
                             Richness_GBIF = rich_gbif)


### Correlation between polygons vs. points of occurence ----
diversity_islands <- left_join(diversity_polygons,
                               diversity_gbif,
                               by = "id")

diversity_islands_cor <- diversity_islands %>%
  remove_missing() 

# Functional dispersion
cor(
  log(diversity_islands_cor$FDisp+1),
  log(diversity_islands_cor$FDisp_GBIF+1),
  method = "pearson"
)
# Phylogenetic endemism
cor(
  log(diversity_islands_cor$PhyEnd),
  log(diversity_islands_cor$PhyEnd_GBIF),
  method = "pearson"
) 
# Richness
cor(
  log(diversity_islands_cor$Richness),
  log(diversity_islands_cor$Richness_GBIF),
  method = "pearson"
) 

### Functional diversity imputed vs. non-imputed ----
# obtain functional diversity non-imputed data
# gawdis distance
gaw_anura_nonimputed <- gawdis(
  amphibia_without_imputed,
  w.type="analytic",
  groups = c(1,
             2,2,
             3,3,3,3,
             4,4,4))

# Community matrix
polygons_community_nonimputed <- polygons_adj %>%
  column_to_rownames("id") %>%
  select(row.names(amphibia_without_imputed)) %>%
  select(sort(names(.))) %>%
  as.matrix()

ncol(polygons_community_nonimputed) # 991 especies
nrow(polygons_community_nonimputed) # 1227 comunidades 

remove <- which(rowSums(polygons_community_nonimputed) == 0)
polygons_community_nonimputed_adj <- polygons_community_nonimputed[-remove, ]

ncol(polygons_community_nonimputed_adj) # 991 especies
nrow(polygons_community_nonimputed_adj) # 1072 comunidades  

# FDis
disp_func_amphibia_nonimputed <- fdisp(
  gaw_anura_nonimputed, 
  polygons_community_nonimputed_adj,
  tol = 1e-07)

df_noninputed <- data.frame(
  id = as.integer(row.names(polygons_community_nonimputed_adj)),
  FDis_nonimputed = disp_func_amphibia_nonimputed$FDis)

df_functional_validation <- diversity_polygons %>%
  select(id, FDisp) %>%
  right_join(df_noninputed, by = "id")

# Save data ----
# Polygons diversity metrics
write.csv2(diversity_polygons, file.path(local_directory, "biogeo_islands",
                                         "amp_islands_diversity.csv"))

# GBIF diversity metrics
write.csv2(diversity_gbif, file.path(local_directory, "biogeo_islands",
                                     "amp_islands_diversity_GBIF.csv"),
           row.names = FALSE)

# Exploring results ----
getwd() # pasta que ira salvar as imagens
ggplot(diversity_polygons, aes(x = FRich, y = log(Richness))) +
  geom_point(size = 3, color = "blue") +  # Pontos azuis
  geom_smooth(method = "lm", se = TRUE, color = "red", fill = "pink") + # Linha de regressão
  labs(x = "FRich", y = "Log(Richness)") +
  theme_minimal() 

diversity_test <- diversity_polygons %>%
  remove_missing()

cor(log(diversity_test$FRich),
    log(diversity_test$Richness), method = "pearson") # 0.75 correlacionadas

ggplot(diversity_test, aes(x = log10(FRich), y = log10(Richness))) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "loess", col = "blue") +
  labs(title = "",
       x = "Functional Richness", y = "Species richness") +
  theme_minimal() 

# Cross-validation
# transform and plot data correlation
df_plot <- diversity_islands_cor %>%
  mutate(
    log_FDisp = log(FDisp + 1),
    log_FDisp_GBIF = log(FDisp_GBIF + 1),
    log_PhyEnd = log(PhyEnd),
    log_PhyEnd_GBIF = log(PhyEnd_GBIF),
    log_Richness = log(Richness),
    log_Richness_GBIF = log(Richness_GBIF)
  )

# Gráfico para Functional Dispersion
p1 <- ggplot(df_plot, aes(x = log_FDisp, y = log_FDisp_GBIF)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", col = "blue") +
  labs(title = "Functional dispersion",
       x = "IUCN", y = "GBIF") +
  theme_minimal()

# Gráfico para Phylogenetic Endemism
p2 <- ggplot(df_plot, aes(x = log_PhyEnd, y = log_PhyEnd_GBIF)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", col = "blue") +
  labs(title = "Phylogenetic endemism",
       x = "IUCN", y = "GBIF") +
  theme_minimal()

# Gráfico para Richness
p3 <- ggplot(df_plot, aes(x = log_Richness, y = log_Richness_GBIF)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", col = "blue") +
  labs(title = "Species Richness",
       x = "IUCN", y = "GBIF") +
  theme_minimal()

tiff("correlation_metrics.tiff", width = 1500, height = 3000, res = 300)
cowplot::plot_grid(p1,p2,p3, ncol = 1, align = "v")
dev.off()

# Validation imputed vs. non imputed data
cor(log(df_functional_validation$FDisp+1),
    log(df_functional_validation$FDis_nonimputed+1),
    method = "pearson") # 0.75 correlacionadas

tiff("correlation_metrics.tiff", width = 1500, height = 3000, res = 300)
ggplot(df_functional_validation, aes(x = log(FDisp+1), 
                                     y = log(FDis_nonimputed+1))) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", col = "blue") +
  labs(title = "Validation inputed and non-imputed traits",
       x = "FDis with inputed traits", y = "FDis with non-inputed traits") +
  theme_minimal() 


tiff("imputed_metrics.tiff", width = 1800, height = 1000, res = 300)
p <- ggplot(df_functional_validation, aes(x = log(FDisp + 1), y = log(FDis_nonimputed + 1))) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", col = "blue") +
  labs(title = "",
       x = "Imputed data", y = "Non-imputed data") +
  theme_minimal()

# Adicionar histogramas marginais
ggMarginal(p, type="density", fill = "gray", alpha = 0.3)
dev.off()