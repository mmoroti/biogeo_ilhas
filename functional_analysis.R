# Mon Aug 22 14:09:37 2022 ------------------------------
# functional analysis for anunan biogegraphy

# functions
# Imputation body_size
# I modified the function fitGeospiza so that it receives the two arguments (phy, trait) that will go in the model. The function's new name is fit_modified2. In fit_modified2, I included the aicw function from the geiger package to also return us the AICw in the model comparison.
#---- CREATE a FUNCTION to COMPARE EVOLUTION MODELS
fit_modified2=function(phy, trait){
  #trait=match.arg(trait, c("body_size","litter_size"))
  
  # define set of models to compare
  models=c("BM", "OU", "EB", "white")
  summaries=c("diffusion", "Ornstein-Uhlenbeck", "early burst", "white noise")
  
  ## ESTIMATING measurement error ##
  aic.se=numeric(length(models))
  lnl.se=numeric(length(models))
  w.se=numeric(length(models))
  for(m in 1:length(models)){
    cat("\n\n\n\n\t*** ", paste(toupper(summaries[m]),": fitting ", sep=""), models[m],
        " with SE *** \n", sep="")
    tmp=fitContinuous(phy, trait,SE=NA, model=models[m],
                      bounds=list(SE=c(0,0.5)), ncores=2)
    print(tmp)
    aic.se[m]=tmp$opt$aicc
    lnl.se[m]=tmp$opt$lnL
    w.se[m]=tmp$opt$aic
  }
  
  
  ## ASSUMING no measurement error ##
  aic=numeric(length(models))
  lnl=numeric(length(models))
  w.se=numeric(length(models))
  
  for(m in 1:length(models)){
    cat("\n\n\n\n\t*** ", paste(toupper(summaries[m]),": fitting ", sep=""), models[m],
        " *** \n", sep="")
    tmp=fitContinuous(phy, trait,SE=0,model=models[m], ncores=2)
    print(tmp)
    aic[m]=tmp$opt$aicc
    lnl[m]=tmp$opt$lnL
    w.se[m]=tmp$opt$aic
  }
  
  ## COMPARE AIC ##
  names(aic.se)<-names(lnl.se)<-names(aic)<-names(lnl)<-models
  delta_aic<-function(x) x-x[which(x==min(x))]
  
  # no measurement error
  daic=delta_aic(aic)
  cat("\n\n\n\t\t\t\t*** MODEL COMPARISON: "," *** \n",sep="")
  cat("\tdelta-AIC values for models assuming no measurement error
    \t\t\t\t zero indicates the best model\n\n")
  print(daic, digits=2)
  
  # measurement error
  daic.se=delta_aic(aic.se)
  cat("\n\n\n\n\t\t\t\t*** MODEL COMPARISON: "," ***\n",sep="")
  cat("\t\t   delta-AIC values for models estimating SE
    \t\t\t\t zero indicates the best model\n\n")
  print(daic.se, digits=2)
  cat("\n\n\n")
  
  # wight akaike
  w=aicw(w.se)
  cat("\n\n\n\n\t\t\t\t*** MODEL COMPARISON: "," ***\n",sep="")
  cat("\t\t   w-AIC values for models estimating SE
    \t\t\t\t zero indicates the best model\n\n")
  print(w, digits=2)
  cat("\n\n\n")
  
  res_aicc=rbind(aic, aic.se, daic, daic.se, w$w)
  rownames(res_aicc)=c("AICc","AICc_SE","dAICc", "dAICc_SE", "AICw")
  #print(w)
  return(res_aicc)
}

# packages
library(FD)
library(tidyverse)
library(picante)
library(RRphylo)
library(Rphylopars)
library(phytools)
library(geiger)

# loading data
#Composition data
dir()
composition <- as_tibble(read.csv2("anura_islands_without0.csv"))
# put names in rows
#View(composition)
composition.select <- composition %>% select(-X, -id)
rownames(composition.select) <- as.character(composition$id) 
View(as.data.frame(rowSums(composition.select)))
# phylogeny
phy_anura_islands <- ape::read.tree("anura_islands.tre")
phy_anura_islands #1986

# traits 
dir()
traits_anura <- as_tibble(read.csv2("anura_islands_traits.csv")) %>% select(-X)
visdat::vis_miss(traits_anura)
nrow(traits_anura)

# preparing continuous variable 'body_size_mm' for imputation data. 12,37% are missing data.
body_size <- as.data.frame(traits_anura %>% select(Species, Body_size_mm) %>% rename(species = Species))

body_size$Body_size_mm <- log10(body_size$Body_size_mm)

#rename rownames
body_size_to_input <- body_size %>% select(-species)

rownames(body_size_to_input) <- body_size$Species
nrow(body_size_to_input)

# We need to remove NA's to do the selection for best fit model of trait evolution
body_size_withoutna <- remove_missing(body_size, vars="Body_size_mm") #238 spp removed
View(body_size_withoutna) #1686
# preparing dataset
body_size_without_rownames <- body_size_withoutna %>% select(-species)
rownames(body_size_without_rownames) <- body_size_withoutna$species 
View(body_size_without_rownames)

# Colnames with species ID - This proceding is necessary for prune.sample function
body_amphi_trans <- t(body_size_without_rownames)
#colnames(body_amphi_trans) <- body_amphi_trans[1,] 
View(body_amphi_trans)
ncol(body_amphi_trans)
# prune sample with phylogeny
phy_body <- prune.sample(body_amphi_trans,phy_anura_islands)
phy_body # 1686 tips
nrow(body_size_without_rownames) #1686 spp

# Aqui inicia o processo de sele????o de modelos evolutivos para os atributos, precisamos usar os dados sem NA's depois fazemos a imputa????o. 
# Calculing the best fit evolution model
# Function fit_modified2 is required
###--- Anura
# The Phylogeny needs to rooted and ultrametric
is.ultrametric(phy_body)
anura_phy_ultra <- force.ultrametric(phy_body)
anura_phy_ultra <- ape::multi2di(anura_phy_ultra)

# Body models
body_anura_model <- fit_modified2(anura_phy_ultra, body_size_without_rownames)

body_anura_model # OU is the best model

#--- after selection evolution model
anura_phy_ultra_full <- force.ultrametric(phy_anura_islands)
amphibia_phy_rooted <- ape::multi2di(anura_phy_ultra_full)
is.ultrametric(amphibia_phy_rooted)

p_OU <- phylopars(trait_data = body_size, tree = amphibia_phy_rooted,  pheno_error = TRUE,phylo_correlated = TRUE,pheno_correlated = TRUE, model="mvOU")
View(p_OU$anc_recon)
phy_bodysize <- as.data.frame(p_OU$anc_recon[1:1924,]) # Data with imputed species means

# Join phylogenetic imputation in datase
body_input <- as_tibble(cbind(body_size_input = phy_bodysize$`p_OU$anc_recon[1:1924, ]`, Species = rownames(phy_bodysize)))
body_input$body_size_input <- as.numeric(body_input$body_size_input)

# Join without NA's in body size
traits_input <- left_join(body_input, traits_anura, by="Species")  %>% select(-Body_size_mm)
hist(traits_input$body_size_input)#normal distribution

# Thu Aug 25 18:22:18 2022 ------------------------------
# traits with new composition 
# Arrumando para linhas
anuran_list <- composition.select %>%
  gather(key="Species",value="count", Abavorana_luctuosa:Zhangixalus_viridis) %>% distinct(Species, .keep_all = TRUE)
nrow(anuran_list)

traits_input <- left_join(anuran_list, traits_input, by="Species")  %>% select(-count)
View(traits_input)

# removing missing data in development mode and habitat
visdat::vis_miss(traits_input)
View(traits_input)
# We only 4.4% percent missing data
# 3.53% habitat (fos,ter,aqu,arb)
# 8.58% development

traits_without_na <- remove_missing(traits_input, vars=names(traits_input))
nrow(traits_without_na) # 1736 species with functional data
#we missing 188 spp

# Wed Aug 24 15:10:03 2022 ------------------------------
# functional metrics 
names(traits_without_na)
# development 
development <- traits_without_na %>% 
  select(Species, Dir, Lar, Viv) %>%
  mutate_at(vars(2:4), list(~replace(., is.na(.), 0)))

development <- column_to_rownames(as.data.frame(development),var = "Species")

dist_development <- prep.binary(development, col.blocks = 3, labels = "development")

# habitat
habitat <- traits_without_na %>% 
  select(Species, Fos, Ter, Aqu, Arb) %>%
  mutate_at(vars(2:5), list(~replace(., is.na(.), 0)))
habitat

habitat <- column_to_rownames(as.data.frame(habitat), 
                              var = "Species")

dist_habitat <- prep.binary(habitat, col.blocks = 4, labels = "habitat")

# Transformando para numeric traitss
sum(is.na(traits_without_na$body_size_input)) # 0 NA's

#traits.amphibia.phy <- traits.amphibia.phy %>% 
#  mutate_at("Body_size_mm", list(as.numeric))
#class(traits.amphibia.phy$Body_size_mm)

#Explorando body size
ggplot(traits_without_na, aes(body_size_input))+
  geom_density() #normal distribution

bodysize <- data.frame(body_size = traits_without_na$body_size_input)
rownames(bodysize) <- traits_without_na$Species
head(bodysize)

# conferindo
nrow(bodysize)#1736
nrow(habitat)
nrow(development)

#Now, let's finally calculate the pair-wise trait distance matrix using the Gower distance coefficient as implemented in the package `ade4'
trait.dist <- dist.ktab(ktab.list.df(list(dist_habitat,dist_development, bodysize)), type = c("B","D","Q"), scan=TRUE)

?dist.ktab
is.euclid(trait.dist)
#10 = S2 coefficient of GOWER & LEGENDRE for multi-choice traits
#1 = Euclidean for quantitative traits

####
# Wed Aug 24 15:25:33 2022 ------------------------------

# At first, we need to construct our composition matrix only species with functional data available 
# Tirando da lista
list <- anti_join(traits_anura,traits_without_na, by="Species")

list.remove <- list$Species

dim(composition.select)

#Removendo da composi????o (altera????o do nome de . para _)
composition_functional <- composition.select[,!(names(composition.select)%in% list.remove)]
dim(composition_functional) #1222 ilhas

# identificando ilha com zero esp??cies
teste <- cbind(id=composition$id, composition_functional)
View(teste)
rownames(teste) <- composition$id

# check match
View(as.data.frame(rowSums(teste[,-1])))
View(as.data.frame(rowSums(composition_functional))) # one community has 0 spp, because this in this functional analysis we removed more one community (island id 6938).

composition_functional <- composition_functional[-665,]
composition_id <- teste[-665,]
View(composition_id)
# Now, we calculating Functional Diversity metrics
disp.func.amphibia <- fdisp(trait.dist, as.matrix(composition_functional),  tol = 1e-07)
disp.func.amphibia

# Histograma
disp.func.amphibia$FDis
hist(disp.func.amphibia$FDis)

# join data frame with id islands
FDis_metric_amphibia <- data.frame(id=composition_id$id, FDis=disp.func.amphibia$FDis) #save csv file with functional metrics 
View(FDis_metric_amphibia)
# fdis
dir()
write.csv2(FDis_metric_amphibia, "fdis_islands.csv")

# Thu Sep 01 18:04:25 2022 ------------------------------
# Functional richness
traits_teste <- traits_without_na
rownames(traits_teste) <- traits_without_na$Species
nrow(traits_teste)#1736

teste_comp <- composition_id[,-1]
View(teste_comp)

ex3 <- dbFD(trait.dist, teste_comp)

#grouping and save
fric_index <- cbind(id=names(ex3$FRic),FRic=(ex3$FRic))
write.csv2(fric_index, "fric_islands.csv")
