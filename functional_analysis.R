# Mon Aug 22 14:09:37 2022 ------------------------------
# functional analysis for anunun biogegraphy
library(FD)
library(tidyverse)
library(picante)
library(RRphylo)
library(Rphylopars)
library(phytools)

# loading data
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

# Aqui inicia o processo de seleção de modelos evolutivos para os atributos, precisamos usar os dados sem NA's depois fazemos a imputação. 
# Calculing the best fit evolution model
# Function fit_modified2 is required
###--- Anura
# The Phylogeny needs to rooted and ultrametric
anura_phy_ultra_2 <- force.ultrametric(anura_phy_body)
amphibia_phy_rooted_2 <- ape::multi2di(anura_phy_ultra_2)

# Body models
body_anura_model <- fit_modified2(amphibia_phy_rooted_2, svl_amph_sem)
body_anura_model # dAICc OU model - AICw 1.0 OU

#--- after selection evolution model
anura_phy_ultra <- force.ultrametric(anura_phy)
amphibia_phy_rooted <- ape::multi2di(anura_phy_ultra)
is.ultrametric(amphibia_phy_rooted)

p_OU <- phylopars(trait_data = short_amphibia_traits, tree = amphibia_phy_rooted,  pheno_error = TRUE,phylo_correlated = TRUE,pheno_correlated = TRUE, model="mvOU")

amphibia_phy_rooted

p_OU$anc_recon[1:20,] # Data with imputed species means
p_OU$anc_var[1:20,] # Variances for each estimate

p_OU$anc_recon[1:20,] - sqrt(p_BM$anc_var[1:20,])*1.96 # Lower 95% CI
p_OU$anc_recon[1:20,] + sqrt(p_BM$anc_var[1:20,])*1.96 # Upper 95% CI
View(p_OU$anc_recon)

# removing missing data in development mode and habitat


# functional metrics 
# Primeiro, precisamos preparar as variáveis
#Atividade circadiana
atividade <- traits.amphibia.phy %>% 
  select(sp, Diu, Crepu, Noc) %>%
  mutate_at(vars(2:4), list(~replace(., is.na(.), 0)))
atividade <- as.tibble(atividade)

atividade <- column_to_rownames(as.data.frame(atividade),var = "sp")

dist_ativ <- prep.binary(atividade, col.blocks = 3, labels = "atividade_circadiana")

# Habitat
names(traits.amphibia.phy)
habitat <- traits.amphibia.phy %>% 
  select(sp, Fos, Ter, Aqu, Arb) %>%
  mutate_at(vars(2:5), list(~replace(., is.na(.), 0)))
habitat

habitat <- as.tibble(habitat)

habitat <- column_to_rownames(as.data.frame(habitat), 
                              var = "sp")
dist_habitat <- prep.binary(habitat, col.blocks = 4, labels = "habitat")

# Transformando para numeric traitss
sum(is.na(traits.amphibia.phy$Body_size_mm)) # 0 NA's
traits.amphibia.phy <- traits.amphibia.phy %>% 
  mutate_at("Body_size_mm", list(as.numeric))
class(traits.amphibia.phy$Body_size_mm)

#Explorando body size
ggplot(traits.amphibia.phy, aes(Body_size_mm))+
  geom_density()

bodysize <- data.frame(BodySize=log(traits.amphibia.phy$Body_size_mm))
rownames(bodysize) <- traits.amphibia.phy$sp
head(bodysize)
nrow(bodysize)

#Explorando range extension
sum(is.na(traits.amphibia.phy$area_m))#ok

range.ext <- data.frame(range.extension=(as.numeric(as.character(traits.amphibia.phy$area_m))))

#Unindo traits quantitativos em um conjunto de dados
quant.traits.amphi <- cbind(bodysize, range.ext)
sum(is.na(quant.traits.amphi)) #pronto!

#Now, let's finally calculate the pair-wise trait distance matrix using the Gower distance coefficient as implemented in the package `ade4'
trait.dist <- dist.ktab(ktab.list.df(list(dist_ativ, dist_habitat, as.data.frame(quant.traits.amphi))), type = c("B","B","Q"), scan=TRUE)
#10 = S2 coefficient of GOWER & LEGENDRE for multi-choice traits
#1 = Euclidean for quantitative traits

?dist.ktab

is.euclid(trait.dist)#TRUE =D
sum(as.matrix(is.na(trait.dist)))#no NAs :)

####
# Wed Mar 31 16:11:06 2021 ------------------------------
# Calculating Functional Diversity metrics
disp.func.amphibia <- fdisp(trait.dist, as.matrix(list.amphibia.phy))

FDis_metric_amphibia <- data.frame(id=grid_ma$id, FDis=disp.func.amphibia$FDis)

# save csv file with functional metrics 