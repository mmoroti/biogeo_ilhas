# Mon Jul  8 12:17:12 2024 ------------------------------
library(phyloregion)
library(phytools)
library(tidyverse)
library(picante)
library(reshape2)
library(FD)

getwd()
dir()

# Cross-validation using points-of-occurence
# list
list_species_gbif <- read_delim("cross_validation/data_gbif.csv",
                                delim = ",")[,-1]

list_species_iucn <- read_delim("cross_validation/iucn_species.csv",
                                delim = ",")[,-1]


# remove duplicates
list_species_gbif <- list_species_gbif %>% 
  distinct() %>%
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

list_species_gbif$species <- list_species_gbif$species %>%
  str_replace_all("Ranoidea", "Litoria") 

list_species_gbif <- list_species_gbif %>%
  mutate(species = recode(species, 
                          "Amnirana_nicobariensis" = "Bijurana_nicobariensis",
                          "Aphantophryne_parkeri" = "Oreophryne_parkeri",
                          "Leptomantis_harrissoni" = "Rhacophorus_harrissoni",
                          "Leptomantis_angulirostris" = "Rhacophorus_angulirostris",
                          "Leptomantis_belalongensis" = "Rhacophorus_belalongensis",
                          "Leptomantis_bimaculatus" = "Rhacophorus_bimaculatus",
                          "Leptomantis_fasciatus" = "Rhacophorus_fasciatus",
                          "Leptomantis_gauni" = "Rhacophorus_gauni",
                          "Leptomantis_rufipes" = "Rhacophorus_rufipes"
  )) 

list_sp_gbif <- list_species_gbif %>%
  select(species) %>%
  distinct(species)

list_sp_iucn <- list_species_iucn %>%
  select(binomial, OBJECTID_1) %>%
  distinct(binomial)

# Trait match
anura_traits <- as_tibble(
  openxlsx::read.xlsx("reboucas_JBI_v2/anura_traits_raoni_2.xlsx")) %>% 
  select(-X1,-X11,-X12)

# GBIF
join_gbif <- left_join(
  list_sp_gbif,
  anura_traits, 
  by = c("species" = "Species")
) %>%
  arrange(species)
#View(join_gbif)

#View(join_gbif %>% rowwise() %>%
#       filter(all(is.na(across(-species)))) %>%
#       ungroup())
#visdat::vis_miss(join_gbif)

# IUCN
join_iucn <- left_join(
  list_sp_iucn,
  anura_traits,
  by = c("binomial" = "Species")
) %>% arrange(binomial)

join_gbif <- remove_missing(join_gbif, vars=names(join_gbif))
join_iucn <- remove_missing(join_iucn, vars=names(join_iucn))

#View(join_iucn %>% rowwise() %>%
#  filter(all(is.na(across(-binomial)))) %>%
#  ungroup())
#visdat::vis_miss(join_gbif)

# verificando O que saiu
anti_join_gbif <- anti_join(
  list_sp_gbif,
  anura_traits, 
  by = c("species" = "Species")
) %>%
  arrange(species) %>%
  distinct(species)

anti_join_iucn <- anti_join(
  list_sp_iucn,
  anura_traits,
  by = c("binomial" = "Species")
)
# gerar lista de especies que sairam
#getwd()
#write.table(anti_join , "list_unmatch_trait_phy.txt") # 281 especies sem traits

# Composition to calculing functional diversity 
# eh necessario que as especies ocorram em pelo menos uma comunidade
# e que as comunidades nao sejam riqueza = 0
# checando
row_sums <- colSums(matrix_iucn)
linhas_com_soma_zero <- matrix_iucn[, row_sums == 0 ]
colnames(linhas_com_soma_zero)

# GBIF
matrix_species_gbif <- as.data.frame(dcast(
  list_species_gbif, OBJECTID_1 ~ species, length)) %>%
  filter(OBJECTID_1 != "17831",
         OBJECTID_1 != "2920",
         OBJECTID_1 != "2940",
         OBJECTID_1 != "2959",
         OBJECTID_1 != "3278",
         OBJECTID_1 != "3372",
         OBJECTID_1 != "17831",
         OBJECTID_1 != "5017") # remove island with 0 spp.
matrix_gbif <- matrix_species_gbif %>% select(join_gbif$species)
rownames(matrix_gbif) <- matrix_species_gbif$OBJECTID_1

# IUCN
list_sp_iucn <- list_species_iucn %>%
  select(binomial, OBJECTID_1)

matrix_species_iucn <- as.data.frame(dcast(
  list_sp_iucn, OBJECTID_1 ~ binomial, length)) %>%
  filter(OBJECTID_1 != "17831",
         OBJECTID_1 != "2920",
         OBJECTID_1 != "2940",
         OBJECTID_1 != "2959",
         OBJECTID_1 != "3278",
         OBJECTID_1 != "3372",
         OBJECTID_1 != "17831",
         OBJECTID_1 != "5017") # remove island with 0 spp.
matrix_iucn <- matrix_species_iucn %>% select(join_iucn$binomial)
rownames(matrix_iucn) <- matrix_species_iucn$OBJECTID_1

# spp que nao estao presentes em nenhuma comunidade da iucn
matrix_iucn <- matrix_iucn %>% select(- "Alsodes_coppingeri", 
                       - "Alsodes_kaweshkari",
                       - "Chaltenobatrachus_grandisonae")

join_iucn <- join_iucn %>% filter(binomial != "Alsodes_coppingeri",
                     binomial != "Alsodes_kaweshkari",
                     binomial != "Chaltenobatrachus_grandisonae")

# Functional diversity
names(join_gbif)
names(join_iucn)

# development 
development_gbif <- join_gbif %>% 
  select(species, Dir, Lar, Viv) %>%
  mutate_at(vars(2:4), list(~replace(., is.na(.), 0))) %>%
  arrange(species)
development_gbif <- column_to_rownames(as.data.frame(development_gbif),var = "species")
development_gbif <- prep.binary(development_gbif, col.blocks = 3, labels = "development")

development_iucn <- join_iucn %>% 
  select(binomial, Dir, Lar, Viv) %>%
  mutate_at(vars(2:4), list(~replace(., is.na(.), 0))) %>%
  arrange(binomial)
development_iucn <- column_to_rownames(as.data.frame(development_iucn),var = "binomial")
development_iucn <- prep.binary(development_iucn, col.blocks = 3, labels = "development")

# habitat
habitat_gbif <- join_gbif %>% 
  select(species, Fos, Ter, Aqu, Arb) %>%
  mutate_at(vars(2:5), list(~replace(., is.na(.), 0))) %>%
  arrange(species)

habitat_gbif <- column_to_rownames(as.data.frame(habitat_gbif), 
                              var = "species")
habitat_gbif <- prep.binary(habitat_gbif, col.blocks = 4, labels = "habitat")

habitat_iucn <- join_iucn %>% 
  select(binomial, Fos, Ter, Aqu, Arb) %>%
  mutate_at(vars(2:5), list(~replace(., is.na(.), 0))) %>%
  arrange(binomial)

habitat_iucn <- column_to_rownames(as.data.frame(habitat_iucn), 
                                   var = "binomial")
habitat_iucn <- prep.binary(habitat_iucn, col.blocks = 4, labels = "habitat")

# body size
bodysize_gbif <- data.frame(body_size = join_gbif$Body_size_mm)
rownames(bodysize_gbif) <- join_gbif$species

bodysize_iucn <- data.frame(body_size = join_iucn$Body_size_mm)
rownames(bodysize_iucn) <- join_iucn$binomial

# conferindo
nrow(bodysize_gbif)
nrow(habitat_gbif)
nrow(development_gbif)

nrow(bodysize_iucn)
nrow(habitat_iucn)
nrow(development_iucn)

#Now, let's finally calculate the pair-wise trait distance matrix using the Gower distance coefficient as implemented in the package `ade4'
trait.dist <- dist.ktab(ktab.list.df(list(habitat_gbif,
                                          development_gbif,
                                          bodysize_gbif)), type = c("B","D","Q"),
                        scan=TRUE)
trait.dist.iucn <- dist.ktab(ktab.list.df(list(habitat_iucn,
                                          development_iucn,
                                          bodysize_iucn)), type = c("B","D","Q"),
                        scan=TRUE)
#10 = S2 coefficient of GOWER & LEGENDRE for multi-choice traits
#1 = Euclidean for quantitative traits

# Now, we calculating Functional Diversity metrics
disp.func.amphibia <- fdisp(trait.dist, as.matrix(matrix_gbif),  tol = 1e-07)

disp.func.amphibia.iucn <- fdisp(trait.dist.iucn,
                                 as.matrix(matrix_iucn),  tol = 1e-07)

FDis_gbif <- data.frame(id=matrix_species_gbif$OBJECTID_1,
                        FDis_gbif =disp.func.amphibia$FDis)

FDis_iucn <- data.frame(id=matrix_species_iucn$OBJECTID_1,
                        FDis_iucn =disp.func.amphibia.iucn$FDis)

join_fdis <- inner_join(FDis_gbif,
                        FDis_iucn,
                        by = "id")

# plot phylogenetic diversity 
ggplot(join_fdis, aes(x = log(FDis_gbif), y = log(FDis_iucn))) +
  geom_point(color = "blue") +
  geom_smooth(method = "lm", col = "red") +
  labs(title = "Functional diversity (FDis)", x = "GBIF", y = "IUCN") +
  theme_minimal()

cor(join_fdis$FDis_gbif, join_fdis$FDis_iucn, method = "pearson")

# Verificar se ha algum zero em cada linha (exceto na primeira coluna)
linhas_sem_zero <- apply(join_fdis[, -1], 1, function(row) all(row != 0))
# Filtrar as linhas que nao contem nenhum valor zero
df_sem_zero <- join_fdis[linhas_sem_zero, ]

ggplot(df_sem_zero, aes(x = log(FDis_gbif), y = log(FDis_iucn))) +
  geom_point(color = "blue") +
  geom_smooth(method = "lm", col = "red") +
  labs(title = "Functional diversity (FDis)", x = "GBIF", y = "IUCN") +
  theme_minimal()
cor(df_sem_zero$FDis_gbif, df_sem_zero$FDis_iucn, method = "pearson")

# Phylogenetic diversity ----
# Takes a community data table and a (rooted) phylogenetic tree (with branch lengths) 
# and calculates either strict or weighted endemism in Phylogenetic Diversity (PD). 
# Strict endemism equates to the total amount of branch length found only in the sample/s
# and is described by Faith et al. (2004) as PD-endemism. Weighted endemism calculates
# the "spatial uniqueness" of each branch in the tree by taking the reciprocal of its range, 
# multiplying by branch length and summing for all branch lengths present at a sample/site.
# Range is calculated simply as the total number of samples/sites at which the branch is present.
# This latter approach is described by Rosauer et al. (2009) as Phylogenetic endemism.
# presence ausence matrix
#matrix_species_gbif <- read_delim("cross_validation/community matrix_GBIF.csv",
#                                  delim = ",")[,-1]
matrix_species_gbif <- as.data.frame(dcast(
  list_species_gbif, OBJECTID_1 ~ species, length))
ncol(matrix_species_gbif)

matrix_species_iucn <- as.data.frame(dcast(
  list_species_iucn, OBJECTID_1 ~ binomial, length))

ncol(matrix_species_iucn)

# Phy match
phy_anura_islands <- ape::read.tree("reboucas_JBI_v2/anura_islands_phy.tre")

match.phylo.comm(phy_anura_islands,matrix_species_gbif)
phy_filter_gbif <- prune.sample(matrix_species_gbif, phy_anura_islands) # 289 spp perdidas
phy_filter_iucn <- prune.sample(matrix_species_iucn, phy_anura_islands)


# filter dataframe
rownames(matrix_species_gbif) <- matrix_species_gbif$OBJECTID_1
rownames(matrix_species_iucn) <- matrix_species_iucn$OBJECTID_1

matrix_gbif <- matrix_species_gbif %>% select(phy_filter_gbif$tip.label)
matrix_iucn <- matrix_species_iucn %>% select(phy_filter_iucn$tip.label)

# alguma comunidade com zero?
any(rowSums(matrix_gbif) == 0)
any(rowSums(matrix_iucn) == 0)

# Analysis of phylogenetic endemism and phyloregions 
#Sparse community matrix - requirement from analysis of phyloregion package
matrix_sparse_gbif <- dense2sparse(matrix_gbif)
matrix_sparse_iucn <- dense2sparse(matrix_iucn)

# phylo_endemism
pe_gbif <- phylo_endemism(matrix_sparse_gbif, phy_filter_gbif)
pe_iucn <- phylo_endemism(matrix_sparse_iucn, phy_filter_iucn)

# 
pe_gbif <- data.frame(ID = rownames(data.frame(pe_gbif)), 
                      pe_gbif)
pe_iucn <- data.frame(ID = rownames(data.frame(pe_iucn)),
                      pe_iucn)

join <- left_join(pe_gbif,
               pe_iucn,
               by="ID")
head(join)

# plot phylogenetic diversity 
ggplot(join, aes(x = log(pe_gbif), y = log(pe_iucn))) +
  geom_point(color = "blue") +
  geom_smooth(method = "lm", col = "red") +
  labs(title = "Phylogenetic endemism (PE)", x = "GBIF", y = "IUCN") +
  theme_minimal()

cor(join$pe_gbif, join$pe_iucn, method = "pearson")

# Richness ----
richness_gbif <- data.frame(ID = matrix_species_gbif[,1],
                            rich_gbif = rowSums(matrix_species_gbif[,-1]))

richness_iucn <- data.frame(ID = matrix_species_iucn[,1],
                            rich_iucn = rowSums(matrix_species_iucn[,-1]))

join_rich <- left_join(richness_gbif,
                       richness_iucn,
                       by="ID")
head(join_rich)

ggplot(join_rich, aes(x = log(rich_gbif), y = log(rich_iucn))) +
  geom_point(color = "blue") +
  geom_smooth(method = "lm", col = "red") +
  labs(title = "Species richness", x = "GBIF", y = "IUCN") +
  theme_minimal()

cor(join_rich$rich_gbif, join_rich$rich_iucn, method = "pearson")
