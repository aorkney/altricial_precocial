# The purpose of this script is to explore the likely evolutionary
# history of divergence between wing and leg skeletal proportions in birds,
# and to determine whether the rate of divergent wing and leg evolution
# differs substantially between bird lineages with high parental investment (altricial)
# and lineages with low parental investment (precocial). 
# We hypothesise that altricial development will alleviate physiological 
# hurdles to postembryonic wing skeleton growth, enabling greater divergence
# and evolutionary rates of change in altricial lineages. 
# We will also explore potential relationships divergence has with 
# body mass variety across birds, as this may indicate an ontogenetic/growth
# mediated mechanism.
# We will also explore whether divergent wing and leg skeletal traits
# are associated with unusual flight-style variety across birds. 


# The McTavish et al. 2025 avian super tree will be used to represent
# evolutionary relationships between birds:
# https://doi.org/10.1073/pnas.2409658122

library(clootl) # 0.1.1
# Load tree package

Tavish.Tree <- extractTree(
species = "all_species",
label_type = "scientific",
taxonomy_year = 2023,
version = "1.4",
data_path = FALSE
)
# Load McTavish super tree

# The curious user may seek out and use their own phylogeny. 
# Popular alternatives might include Prum et al., 2015, Jetz et al., 2012 or Stiller et al., 2024.
# We caution the user that these phylogenies will not include all bird taxa. 
# They may either need to map bird species we investigate to appropriate tips, and investigate smaller datasets, 
# or use an algorithm to replace the backbone of the McTavish tree with the backbones available from these other topologies. 

setwd('') # Ammend to your filepath where you downloaded the data. 

Bjarnason <- read.csv('Ontogeny_Bjarnason_04_30_2025.csv')
# Read Bjarnason's bone length data, the ontogenetic categorizations and broader clade affiliations.
# I sourced the original measurements from Bjarnason et al., 2021: https://doi.org/10.18563/journal.m3.125
# Developmental categorizations are sourced from Starck et al., 1993 https://doi.org/10.1007/978-1-4615-9582-3_6

Bj.birds <- paste(Bjarnason$Genus,Bjarnason$Species,sep='_')
# Original museum names associated with the specimens. 
# Some names have become deprecated over time and we need to update them. 

#'Anas_discors' should be  'Spatula_discors'
#'Apaloderma_narino' should be 'Apaloderma_narina' 
#'Choriotis_kori' should be 'Ardeotis_kori' %in% Tavish.Tree$tip
#'Caracara_cheriway' should be 'Caracara_plancus' 
#'Pipra_erythrocephala' should be 'Ceratopipra_erythrocephala' %in% Tavish.Tree$tip
#'Charadris_vociferus' should be 'Charadrius_vociferus' %in% Tavish.Tree$tip
#'Choloroceryle_amazona' should be 'Chloroceryle_amazona'%in% Tavish.Tree$tip
#'Larus_novaehollandiae' should be 'Chroicocephalus_novaehollandiae' %in% Tavish.Tree$tip
#'Climacteris_melanura' should be 'Climacteris_melanurus' %in% Tavish.Tree$tip
#'Columbina_minutas' should be 'Columbina_minuta' %in% Tavish.Tree$tip
#'Crax_mitu' should be 'Mitu_mitu' 
#'Chloropipio_holochlora' should be 'Cryptopipo_holochlora'
#'Cuculus_fugax' should be 'Hierococcyx_fugax'
#'Diomedea_irrorata' should be 'Phoebastria_irrorata'
#'Elanus_caerulenus' should be 'Elanus_caeruleus'
#'Grus_leucogeranus' should be 'Leucogeranus_leucogeranus'
#'Hirundinea_bellicosa' should be 'Hirundinea_ferruginea'
#'Hymenops_persipicillata' should be 'Hymenops perspicillatus'
#'Leptoptilos_crumeniferus' should be 'Leptoptilos_crumenifer'
#'Megalaima_chrysopogon' should be 'Psilopogon_chrysopogon'
#'Oceanodroma_leucorhoa' should be 'Hydrobates_leucorhous'
#'Odontophorous_guttatus' should be 'Odontophorus_guttatus'
#'Pelecanoides_urinatrix' should be 'Pelecanoides_urinatrix'
#'Phalacrocorax_albiventer' should be 'Leucocarbo_atriceps'
#'Ptilinopus_lechlancheri' should be 'Ptilinopus_leclancheri'
#'Puffinus_tenuirostris' should be 'Ardenna_tenuirostris'
#'Pycnonotus_caffer' should be 'Pycnonotus_cafer'
#'Rallus_limnicola' should be 'Rallus_limicola'
#'Rhamphastos_ambiguus' should be 'Ramphastos_ambiguus'
#'Regulus_ignicapillus' should be 'Regulus_ignicapilla'
#'Rollandia_rollandia' should be 'Rollandia_rolland'
#'Terenura_spodioptila' should be 'Euchrepomis_spodioptila'
#'Tockus_nasutus' should be 'Lophoceros_nasutus'
#'Tringa_ochrophus' should be 'Tringa_ochropus'

Bjarnason$Tavish <- Bj.birds
Bjarnason$Tavish[which(Bj.birds%in%Tavish.Tree$tip==F)] <- c('Spatula_discors',
'Apaloderma_narina','Ardeotis_kori','Caracara_plancus','Ceratopipra_erythrocephala',
'Charadrius_vociferus','Chloroceryle_amazona','Chroicocephalus_novaehollandiae',
'Climacteris_melanurus','Columbina_minuta','Mitu_mitu','Cryptopipo_holochlora',
'Hierococcyx_fugax','Phoebastria_irrorata','Elanus_caeruleus','Leucogeranus_leucogeranus',
'Hirundinea_ferruginea','Hymenops_perspicillatus','Leptoptilos_crumenifer',
'Psilopogon_chrysopogon','Hydrobates_leucorhous','Odontophorus_guttatus',
'Pelecanoides_urinatrix','Leucocarbo_atriceps','Ptilinopus_leclancheri',
'Ardenna_tenuirostris','Pycnonotus_cafer','Rallus_limicola','Ramphastos_ambiguus',
'Regulus_ignicapilla','Rollandia_rolland','Euchrepomis_spodioptila','Lophoceros_nasutus','Tringa_ochropus')

# Name updated complete. 

short.brinkworth <- read.csv('ontogeny_magnificent.csv')
# Read another dataset of bird bone length measurements, developmental strategies etc.
# I sourced and quality controlled the dataset from: https://doi.org/10.1038/s41467-023-41415-2
# Developmental categorizations are sourced from Starck et al., 1993 https://doi.org/10.1007/978-1-4615-9582-3_6

# As before, deprecated species names must be updated to contemporary names in the McTavish super tree. 

short.brinkworth$Tavish <- short.brinkworth$Species
short.brinkworth$Tavish[which(short.brinkworth$Tavish%in%Tavish.Tree$tip==F)] <- c('Gymnogyps_californianus',
'Aechmophorus_occidentalis','Duplicate','Nannopterum_auritum','Leucocarbo_bougainvillii', 'Urile_urile',
'Botaurus_stellaris', 'Platalea_leucorodia','Lophonetta_specularioides','Uncertain_species', 'Daptrius_australis',
'Gypaetus_barbatus','Uncertain_species','Milvus_migrans','Esacus_magnirostris','Stercorarius_skua',
'Leucophaeus_atricilla','Larus_glaucescens', 'Duplicate', 'Saudareos_ornata', 'Callipepla_gambelii',
'Duplicate', 'Strix_virgata','Todiramphus_sanctus','Buceros_rhinoceros','Duplicate','Pitangus_sulphuratus',
'Dicrurus_macrocercus','Corcorax_melanorhamphos','Mareca_americana','Spatula_cyanoptera','Duplicate',
'Mareca_falcata','Spatula_hottentota','Mareca_sibilatrix','Spatula_smithii','Mareca_strepera',
'Spatula_versicolor','Eupetomena_cirrochloris','Eupsittula_canicularis','Psittacara_erythrogenys',
'Psittacara_mitratus','Eupsittula_nana','Geranoaetus_polyosoma','Eolophus_roseicapilla',
'Stercorarius_maccormicki','Chalcopsitta_scintillata','Buteogallus_solitarius','Leucophaeus_fuliginosus',
'Leucophaeus_modestus','Duplicate','Chroicocephalus_philadelphia','Chroicocephalus_serranus',
'Larus_glaucoides','Neophema_bourkii','Nannopterum_brasilianum','Flightless','Urile_pelagicus',
'Urile_penicillatus','Glossoptilus_goldiei','Ardenna_bulleri','Ardenna_carneipes',
'Ardenna_creatopus','Ardenna_gravis','Ardenna_grisea','Ardenna_pacifica', 'Duplicate',
'Phalaropus_tricolor','Onychoprion_anaethetus','Sternula_antillarum','Thalasseus_bergii',
'Hydroprogne_caspia','Thalasseus_elegans','Onychoprion_lunatus','Radjah_radjah', 'Duplicate',
'Anarhynchus_montanus','Chroicocephalus_ridibundus','Pampusana_beccarii','Spilopelia_chinensis',
'Spilopelia_senegalensis','Anser_rossii','Oressochen_melanopterus',
'Pternistis_adspersus','Pternistis_afer','Pternistis_leucoscepus','Oressochen_jubatus',
'Zapornia_flavirostra','Pterorhinus_albogularis','Petrochelidon_spilodera',
'Daptrius_chimachima','Psephotellus_varius','Calidris_virgata','Tringa_semipalmata','Tringa_incana',
'Aramides_cajaneus')

# Updates done



short.brinkworth$Species <- short.brinkworth$Tavish
short.brinkworth <- short.brinkworth[,1:21]
# Book keeping
 
combined <- short.brinkworth
for(i in 1:length(Bjarnason$Tavish)){
	if(Bjarnason$Tavish[i] %in% short.brinkworth$Species){
		combined[match(Bjarnason$Tavish[i], short.brinkworth$Species),]<- t(as.matrix(c(short.brinkworth[match(Bjarnason$Tavish[i], short.brinkworth$Species),1],Bjarnason[i,]$Order,Bjarnason[i,]$Family,Bjarnason$Tavish[i],short.brinkworth[match(Bjarnason$Tavish[i], short.brinkworth$Species),]$Jetz,Bjarnason[i,]$humerus,
			NA,Bjarnason[i,]$radius, Bjarnason[i,]$carpometacarpus, Bjarnason[i,]$femur, Bjarnason[i,]$tibiotarsus, Bjarnason[i,]$tarsometatarsus, NA, 
			short.brinkworth[match(Bjarnason$Tavish[i], short.brinkworth$Species),14], mean(Bjarnason[i,]$Mass_M_.hbw_alive.,Bjarnason[i,]$Mass_F_.hbw_alive.),short.brinkworth[match(Bjarnason$Tavish[i], short.brinkworth$Species),16], short.brinkworth[match(Bjarnason$Tavish[i], short.brinkworth$Species),17],short.brinkworth[match(Bjarnason$Tavish[i], short.brinkworth$Species),18],short.brinkworth[match(Bjarnason$Tavish[i], short.brinkworth$Species),19],Bjarnason[i,]$ontogeny,short.brinkworth[match(Bjarnason$Tavish[i], short.brinkworth$Species),21]   )))
	} else {
			new.line <- t(as.matrix(c(NA,Bjarnason[i,]$Order,Bjarnason[i,]$Family,Bjarnason$Tavish[i],Bjarnason$Tavish[i],Bjarnason[i,]$humerus,
			NA,Bjarnason[i,]$radius, Bjarnason[i,]$carpometacarpus, Bjarnason[i,]$femur, Bjarnason[i,]$tibiotarsus, Bjarnason[i,]$tarsometatarsus, NA, 
			NA, mean(Bjarnason[i,]$Mass_M_.hbw_alive.,Bjarnason[i,]$Mass_F_.hbw_alive.), NA, NA,NA,NA,Bjarnason[i,]$ontogeny,NA   )))
			colnames(new.line) <- colnames(combined)
			combined <- rbind(combined,new.line)
	}
}

# These loops combine the two datasets from Bjarnason and Brinkworth. Where a species is repeated, we favour
# Bjarnason's measurements. 
# This results in a sizeable dataset of about 7% of all recognised bird species; hundreds and hundreds. 

# ape_5.8-1 must be installed. 

short.tree <- ape::keep.tip(Tavish.Tree, combined$Species[which(combined$Species%in%Tavish.Tree$tip==T)])
short.tree <- ape::drop.tip(short.tree, 'Acrocephalus_aequinoctialis')
# Prune the grafted tree to the taxa of interest

C <- ape::vcv.phylo(short.tree)
# Shared ancestry matrix
# Brinkworth's data has unusually strong differences in traits between species that diverged under 2ma, suggesting that this horizon
# represents a point where measurement error confounds macroevolutionary differences between species
# We are therefore going to prune the tree so that we retain only a single species within each affected genus:

total.time <- max(C)
C[which(C==diag(C))] <- NA
close.Brinkworth <- names(which(apply(C,1,max,na.rm=T)-total.time >-2))
crop.Brinkworth <- close.Brinkworth 
genera <- unique(gsub('_.*','',crop.Brinkworth))
keeps <- list()
for(i in 1:length(genera)){
	if(genera[i] %in% Bjarnason$Tavish==F){
		keeps[[i]]<-grep(genera[i],crop.Brinkworth)[1]
	}
}
keeps <- unlist(keeps)
crop.Brinkworth <- crop.Brinkworth[-keeps]
crop.Brinkworth <-crop.Brinkworth[which(crop.Brinkworth %in% Bjarnason$Tavish==F)]
# We don't want to lose out on any Bjarnason birds, because we have detailed ecological metadata, so we're going to retain these
# species preferentially. 

# We are cropping birds that have diverged within the last 2 million years, because several bird clades
# experience hybridization between closely related species, especially ducks and geese, which could 
# result in animals with hybrid traits distorting assessments of evolutionary rates.
# Measurement precision may also become a problem at timescales that become very short. 


# Our study concerns hypotheses that exist within the context that all birds share a broadly similar baseline anatomical function
# and set of biomechanical demands on their wings and legs that is at least somewhat analogous. 
# We therefore need to omit flightless species:

flightless.Brinkworth <- c(
'Struthio_camelus',
'Dromaius_novaehollandiae',
'Casuarius_casuarius', 'Casuarius_bennetti',
'Apteryx_owenii',
'Rhea_americana',
'Anas_aucklandica', 
'Tachyeres_brachypterus',
'Gallirallus_australis',
'Gallirallus_owstoni',
'Aptenodytes_patagonicus',
'Pygoscelis_adeliae',
'Eudyptula_minor',
'Spheniscus_demersus',
'Spheniscus_magellanicus',
'Spheniscus_mendiculus',
'Spheniscus_humboldti',
'Megadyptes_antipodes',
'Eudyptes_chrysolophus',
'Eudyptes_chrysocome',
'Phalacrocorax_harrisi','Nannopterum_harrisi')
# We want to remove the flightless birds because our hypothesis is centered on the idea that
# flying birds must navigate the complex challenge of maintaining a functional aerofoil while
# also achieving adaptations to new ecologies- and we hypothesize this task is modulated by
# developmental strategy, which may decouple wing and leg evolution in altricial birds
# and afford access to a new axis of evolutionary innovation 
# Flightless birds, by contrast, are players in a fundamentally different game; they can 
# down regulate or even deactivate their wing-forming developmental programmes, so the origin of
# flightlessness is often associated with rapid wing reductions and departures from concordant wing
# and leg evolution. These birds do not tell us whether developmental strategy modulates solutions
# to adaptation within the context of the biomechanical constraints of flight. 

crop.Brinkworth <- c(crop.Brinkworth,flightless.Brinkworth,'Argusianus_argus')
# I am also removing the Argus, because I suspect it was not measured correctly. 

short.tree <-  ape::drop.tip(short.tree, crop.Brinkworth)
# This is a pruned phylogeny 

C <- ape::vcv.phylo(short.tree)
# Shared ancestry matrix

# Before we investigate bird long bone length variety, we wish to remove strong trends associated with 
# allometry (scale dependent) patterns, because these are likely to overwhelm evolutionary associations evident
# between different bones. 

mass <- log10(as.numeric(combined$Body_mass_g[match(short.tree$tip,combined$Species)]))
names(mass) <- short.tree$tip
# Vector of log-10 transformed body mass

ones<-rep(1,dim(C)[1])
x <- cbind(ones,mass*ones)
# The design matrix for least squares regression of body measurements upon mass 

y=combined[match(short.tree$tip.label,combined$Species),c(6,8,9,10,11,12)]
y <- apply(y,2,as.numeric)
# extract the phenotypic data
y <- apply(y,2,function(x) log10(x) ) 
# Logged limb bone lengths; the dependent dataset
rownames(y) <- short.tree$tip.label

invC <- solve(C); B <- solve(t(x)%*%invC%*%x,tol = 1e-17)%*%t(x)%*%invC%*%(y)
# Find inverse of shared ancestry and solve for allometric trends in limb bone proportions across flying birds

fitted <- t( B[2,]%*%t(x[,2]) + B[1,]%*%t(x[,1]) ); colnames(fitted)<-colnames(y)
# Find the fitted values

r <- y-fitted
# Find residuals from allometric trends



# We are not going to proceed to define 'evolutionary divergence' between wings and legs. 
# We will define this as the orthogonal complement to the vector describing the strongest suite of evolutionary covariances between
# the allometric residuals of wing and leg long bone lengths. 

R<- t(r)%*%invC%*%(r)/dim(C)[1]
# Find covariance matrix among residuals across family tree
pls <- svd(R[1:3,4:6])
# Find major axis of concordant evolution between wing and leg covariance
# c(pls$u[,1],pls$v[,1]) #   -0.5821482 -0.5937930 -0.5554398 -0.5386062 -0.5464676 -0.6413084
# The leading salience of correlated limb bone evolution. Bones share the same sign of coefficient.
# Higher magnitudes indicate variation in this bone's length is more strongly related to correlated 
# evolution between the wing and leg. The Carpometacarpus is not strongly indicative, for example. 
# pls$d/sum(pls$d) # the overwhelming majority of varience is explained by the first axis. (0.980511422)
# Perhaps terms should be squared, makes little difference though

XScores <- r[,1:3] %*% pls$u
YScores <- r[,4:6] %*% pls$v
# Project residuals onto latent vectors describing concordance evolution, to discover
# how bird limb evolution accords with a concordant axis of wing-leg evolution

# Now we are ready to begin thinking about defining sets of birds that differ by developmental strategy: 
development <- rep(NA,dim(combined)[1])
development[grep('tricial_1|tricial_2',combined$Ontogeny)]<-'Altricial'
development[grep('recocial_1|recocial_2|recocial_3|uperrecocial',combined$Ontogeny)]<-'Precocial'
names(development)<-combined$Species
development<-development[short.tree$tip.label]
# Define a vector delineating birds which are strongly precocial or altricial in their developmental style
# Sourced from Starck 1993

reorient <- solve(svd(cbind(XScores[,1],YScores[,1]))$v)
# Rotate the distribution of accorded wing-leg evolution so that the second axis described departure from concordance
# reorient %*% t(cbind(pls$u[,1],pls$v[,1])) # Inspect reoriented coefficients if desired 
#           [,1]       [,2]        [,3]
#[1,]  0.7752074  0.7883293  0.84396839
#[2,] -0.1674712 -0.1724931 -0.08664417
# The reoriented coefficients demonstrate that the concordance axis is dominated by correlated
# evolution between the wing and more proximal components of the leg, and in particular that the 
# zeugopodial component of the wing (radius) is strongly related to the leg. 


divergence <- t(reorient %*% t(cbind(XScores[,1],YScores[,1])))
rownames(divergence) <- short.tree$tip
trait <- as.matrix(divergence[,2])
# Define our axis of divergent wing and leg evolution

short.tree <- phytools::force.ultrametric(short.tree)
# Ensure the tree is ultrametric 

# We will need the following packages: 
library(mvMORPH) # 1.2.1
library( ggplot2 ) # 3.4.1
library( ggdendro ) # 0.1.23
library( dendextend ) # 1.17.1
library(zoo) # 1.8-12
library(dplyr) # 1.1.1

shortest <- drop.tip(short.tree,names(development[which(is.na(development)==T)]))
# Prune our tree to only the developmental sets we are interested in. 

set.seed(3)
# The user is encouraged to change the seed to verify for themselves that different stochastic differences in our code 
# implementation do not affect end results. 

# You will need phytools_2.4-4
simmap.dev <- phytools::make.simmap(tree=shortest, x=development[shortest$tip.label], model="ER", nsim=1,pi=c(0,1))
# Use a stochastic character mapping approach to simulate the possible evolutionary history of developmental strategy across birds;
# Precocial development is assumed as the root state. 


set.seed(3)
# The user is encouraged to change the seed to verify for themselves that different stochastic differences in our code 
# implementation do not affect end results. 

Single <- mvBM(tree=shortest, data= trait[shortest$tip,] , model='BM1')
# Fit a single Brownian rate model of evolutionary divergence between the wing and leg across birds that does not depend on developmental strategy.

set.seed(3)
# The user is encouraged to change the seed to verify for themselves that different stochastic differences in our code 
# implementation do not affect end results. 

Dev.model <- mvBM(tree=simmap.dev, data= trait[shortest$tip,] , model='BMM')
# Model the evolution of divergent wing and leg evolution across birds; is the rate of divergent evolution significantly different
# between altricial and precocial birds? 

LRT(Single,Dev.model)
# Likelihood ratio test 
# P value <<0.001 preferred multiple rates to one rate 
# Inferred rates of 0.000370499 in Altricials verses 0.0001983546 in Precocials; a factor of 1.9 difference to 2SF

# This is an aesthetic function which takes mvMORPH models and produces plots showing the exploration of a single continuous trait
# 'Phenogram' 

make.dendr <- function(tree, trait, mod, obj ){
	dendr <- dendro_data(as.dendrogram(tree))
	lab.dat <- dendr$labels
	fit<-estim(tree, data= trait[tree$tip], object=mod, asr=TRUE)
	dendr.mod<-dendr$segments/2 # This is the tree structure 
	tips <- which(dendr$segments$yend==0)
	nodes <- which(dendr$segments$y==dendr$segments$yend)
	internal.branches <- which(dendr$segments$x==dendr$segments$xend & dendr$segments$yend!=0)
	dendr.mod$x[tips] <- fit$estim[ match( tree$edge[,1][match(dendr$segments$xend[tips],tree$edge[,2])],rownames(fit$estim) )  ] 
	dendr.mod$xend[tips] <- trait[dendr$segments$xend[tips]]
	# This has aligned the tips so that their heads match observation and their tails match ancestral estimates
	dendr.mod[nodes[seq(1,length(nodes),2)] ,]$x <- fit$estim
	dendr.mod[nodes[seq(1,length(nodes),2)]+1 ,]$x <- fit$estim
	dendr.mod$xend[internal.branches] <- dendr.mod$x[internal.branches+3]
	# This gets all of the tips, their nodes and half of the internal branches in place
	dendr.mod[nodes[seq(1,length(nodes),2)]+3 ,]$x <- fit$estim 
	#dendr.mod[nodes[seq(1,length(nodes),2)]+3 ,]$xend <- dendr.mod[nodes[seq(1,length(nodes),2)]+3 ,]$xend
	#dendr.mod$xend[internal.branches] <- fit$estim
	dendr.mod$x[nodes] <- dendr.mod$xend[nodes] <-NA


	for(i in internal.branches){
		# Conditions
		x.match <- which(dendr$segments$x==dendr$segments$xend[i] )
		# which tails of other lines have the same horizontal position as the head of the branch?
		y.match <- which(dendr$segments$y==dendr$segments$yend[i] )
		# which tails of other lines have the same vertical position as the head of the branch?
		head.nodes <- intersect(x.match,y.match)
		# These are the nodes that lead to the next vertical branches
		head.nodes.x.match <- match(dendr$segments$xend[head.nodes], dendr$segments$x)
		head.nodes.x.match <- head.nodes.x.match[which(head.nodes.x.match%in%y.match==T)]
		# which tails of vertical branches meet the head nodes?
		dendr.mod$xend[i] <- dendr.mod$x[head.nodes.x.match][1]	
	}

	dendr.mod$col[tips] <- obj$states[match(tree$edge[,1][match(dendr$segments$xend[tips],tree$edge[,2])],names(obj$states) )]
	dendr.mod[nodes[seq(1,length(nodes),2)] ,]$col <- obj$states
	dendr.mod[nodes[seq(1,length(nodes),2)]+1 ,]$col <- obj$states
	dendr.mod[nodes[seq(1,length(nodes),2)]+2 ,]$col <- obj$states
	dendr.mod[nodes[seq(1,length(nodes),2)]+3 ,]$col <- obj$states
	# Very simply coloring; parent nodal values are inherited. 
	# This works for very simple simmap 

	Rate <- abs((dendr.mod$xend-dendr.mod$x))/(dendr.mod$y-dendr.mod$yend)
	return(cbind(dendr.mod,Rate))
}


obj <- summary(simmap.dev)


phenogram <- make.dendr(tree=simmap.dev, trait=trait[shortest$tip.label,],mod=Dev.model, obj)
# Make the phenogram object 

library(reshape2) # 1.4.4 
# For data formatting and management

dendr <- dendro_data(as.dendrogram(shortest))
lab.dat <- dendr$labels
# Metadata required for plotting 

# Define some plot aesthetics for later use:

cbp <-c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# Define a palette of colours that are friendly to most human visual differences
lwd.factor <- 10
lwd.line <- 1/2



# We now seek to understand whether flight-style variety across birds has an interesting relationship to 
# evolutionary divergence between the wing and leg. 
# We will use the flight-style categorizations compiled by Orkney et al., 2021: https://doi.org/10.1038/s41559-021-01509-w

Flight.style <- read.csv('flight_masses_22_10_2022_plus_A.csv')
# Load ecological classifications for Bjarnason's birds
# We want to determine whether divergent wing and leg evolution is associated with unusual flight style properties
# and whether any such pattern is modulated by developmental strategy

binary.flight.scores <- Flight.style[c(3:13)] ; rownames(binary.flight.scores)<-Flight.style$X
rownames(binary.flight.scores) <- Bjarnason$Tavish[match(rownames(binary.flight.scores),Bjarnason$Prum)]
# Book keeping 

mat <- binary.flight.scores[shortest$tip.label,]
mat <- mat[which(rowSums(mat)>1),]
# Define an ordered matrix for the flight style scores
var<-apply(mat,2,sd)
# Check that columns are not invariant 
mat <- mat[,which(var>0)]
# Remove invariant columns 
flight.dist <- cluster::daisy(mat, metric="gower", type=list('asymm'=1:dim(mat)[2]) )
# transform the flight scores into a distance matrix 
if(length(which(is.na(flight.dist)==T))>0){
	flight.dist[is.na(flight.dist)]<-0
}
# adjust the diagonal to 0 
flight.pcoa <- cmdscale(flight.dist, eig=T)
# Find the explained variances of the eigenvectors of this object
k<-max(which(flight.pcoa$eig/sum(flight.pcoa$eig)>=0.05))
# Consider only those axes which explain more than 5% of variance
flight.pcoa <- cmdscale( flight.dist, k=k  ,eig = T)
# Compute a principal coordinate analysis- a multivariate ordination of the flight style data, constrained to the axes of interest

S<-cov(flight.pcoa$points[,1:k])
center <- colMeans(flight.pcoa$points[,1:k])
mahab <- mahalanobis(flight.pcoa$points[,1:k], center, S)
# Calculate Mahalanobis distance of individual birds' flight style properties from the 'mean' bird.
# We are assuming here that a 'mean' bird is a meaningful concept. This is equivalent to positing that a centerground of 'conventional' flight style
# combinations exists in birds. We are also hypothesizing that birds' flight style evolution tends to emphasize concordant acquisition of flight styles
# in conventional combinations. We are interested in birds which differ strongly from these conventions; do they have unusually divergent limb evolution? 



# Produce a barplot to visualize estimated divergent wing-leg rates in altricial and precocial birds 
bars<- 
	ggplot()+
 	geom_rect(aes(xmin=0,xmax=1,ymin=0,ymax=Dev.model$sigma[,,1:2]['Precocial']),fill='#56B4E9',col='blue',lwd=lwd.line,lineend='round')+
 	geom_rect(aes(xmin=1.2,xmax=2.2,ymin=0,ymax=Dev.model$sigma[,,1:2]['Altricial']),fill='#E69F00',col='yellow',lwd=lwd.line,lineend='round')+
	scale_y_continuous(position = "right", limits=c(0,4.5e-4) )+
	ggpubr::geom_bracket(xmin = 0.5, xmax = 1.7,
    	y.position =  max(Dev.model$sigma[,,1:2])*1.03, label = '***',
    	tip.length = 0.01,size=lwd.line,label.size=11/2,vjust=.55)+
	labs(y=expression(sigma^2))+
	ggtitle("") +
	geom_text(aes(x=c(0.5,1.7),y=c(Dev.model$sigma[,,1:2]['Precocial'],Dev.model$sigma[,,1:2]['Altricial'])/2),label=c('Precocial','Altricial'),angle=90,size=11/.pt)+
	#scale_y_continuous(name = NULL, sec.axis = sec_axis(~., name = expression(sigma^2)))+
	theme( panel.background = element_rect(fill=NA),
         plot.background = element_rect(fill=NA, color=NA),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
		legend.background = element_rect(fill='transparent'),
		legend.key.width=unit( 0.4,'cm'),
		legend.key.height=unit( 0.5,'cm'),
		legend.text=element_text(size=11/2.5,color='black' , family = "Times"),
		axis.line.y.right=element_line( colour = "black",linewidth=lwd.line),
		axis.line.x=element_blank(),
		axis.ticks.x=element_blank(),
		axis.ticks.y.right=element_line(color='black',linewidth=lwd.line),
		axis.title.x=element_blank(),
		axis.title.y.right=element_text(size=40/2.5, color='black',angle=0,vjust=0.95),
		axis.text.x=element_blank(),
		axis.text.y.right=element_text(size=4*.pt,color='black',),
	)


names<-intersect(rownames(flight.pcoa$points),shortest$tip.label)
# Book keeping

flight.dist.df <- data.frame((trait[names,]),mahab[names]^.5,development[names])
colnames(flight.dist.df)<-c('divergence','distinctivenes','development')
linear.altricial <-lm(distinctivenes~divergence + I(divergence^2),data=flight.dist.df[which(flight.dist.df$development=='Altricial'),])
ribbon.altricial.df <-(seq(min(-trait[names(which(development=='Altricial')),]),max(-trait[names(which(development=='Altricial')),]),length.out=80))
ribbon.altricial.df <- data.frame(ribbon.altricial.df) ; colnames(ribbon.altricial.df)<-'trait' ; axis.altricial <- data.frame(ribbon.altricial.df[,1]); colnames(axis.altricial)<-'divergence'
ribbon.altricial.df[,2] <- abs(ribbon.altricial.df[,1])^.5 ; colnames(ribbon.altricial.df)[2]<-'divergence'
ribbon.altricial.df[,3]<-  predict(linear.altricial, newdata = axis.altricial, se.fit = TRUE)$fit
ribbon.altricial.df[,4]<- ribbon.altricial.df[,3]- predict(linear.altricial, newdata = axis.altricial, se.fit = TRUE)$se.fit*1.96
ribbon.altricial.df[,5]<- ribbon.altricial.df[,3]+ predict(linear.altricial, newdata = axis.altricial, se.fit = TRUE)$se.fit*1.96
colnames(ribbon.altricial.df)[3:5]<-c('y','ymin','ymax')
ribbon.altricial.df[,3:5]<- 0.25+(ribbon.altricial.df[,3:5]- min(mahab[rownames(flight.pcoa$points)]^.5))
ribbon.altricial.df[,3:5]<- 35*(ribbon.altricial.df[,3:5]/max(mahab[rownames(flight.pcoa$points)]^.5))
# We hypothesize that the rooted Mahalanobis distance will increase as a quadratic function of divergent wing and leg evolution; this 
# undertakes such a test in altricial birds and produces plotting metadata 
# summary(linear.altricial)
# An OLS fit of rooted Mahalanobis distance of flight style variety as a second order polynomial of divergence wing and leg evolution
# shows that there is a significant dependence on the quadratic term- p <0.01, but that firsr order term does not differ
# signficantly from zero- p ~ 0.83. This is consistent with bird flight style variety exhibitting a 'center ground' of typical ecological
# variety, which is associated with highly correlated wing and leg evolution. By contrast, unusual modes of uncorrelated wing and leg
# evolution are associated with rapidly more divergent flight style combinations. 
# Birds with unusual body plans and evolutionary stories fly in unusual ways, at least within altricial lineages. 


linear.precocial <-lm(distinctivenes~divergence + I(divergence^2),data=flight.dist.df[which(flight.dist.df$development=='Precocial'),])
ribbon.precocial.df <-(seq(min(-trait[names(which(development=='Precocial')),]),max(-trait[names(which(development=='Precocial')),]),length.out=80))
ribbon.precocial.df <- data.frame(ribbon.precocial.df) ; colnames(ribbon.precocial.df)<-'trait' ; axis.precocial <- data.frame(ribbon.precocial.df[,1]); colnames(axis.precocial)<-'divergence'
ribbon.precocial.df[,2] <- abs(ribbon.precocial.df[,1])^.5 ; colnames(ribbon.precocial.df)[2]<-'divergence'
ribbon.precocial.df[,3]<-  predict(linear.precocial, newdata = axis.precocial, se.fit = TRUE)$fit
ribbon.precocial.df[,4]<- ribbon.precocial.df[,3]- predict(linear.precocial, newdata = axis.precocial, se.fit = TRUE)$se.fit*1.96
ribbon.precocial.df[,5]<- ribbon.precocial.df[,3]+ predict(linear.precocial, newdata = axis.precocial, se.fit = TRUE)$se.fit*1.96
colnames(ribbon.precocial.df)[3:5]<-c('y','ymin','ymax')
ribbon.precocial.df[,3:5]<- 0.25+(ribbon.precocial.df[,3:5]- min(mahab[rownames(flight.pcoa$points)]^.5))
ribbon.precocial.df[,3:5]<- 35*(ribbon.precocial.df[,3:5]/max(mahab[rownames(flight.pcoa$points)]^.5))
# We hypothesize that the rooted Mahalanobis distance will not increase as a quadratic function of divergent wing and leg evolution in
# precocial birds, whose wing-leg divergence is limited and whose flight style evolution may be circumscribed within a more conservative space. 
# summary(linear.precocial)
# p ~ 0.73 in precocials for first order and p ~ 0.55 for quadratic, showing that
# uncorrelated modes of wing and leg evolution in precocial birds do not predict unusual flight styles. 

library(cowplot) # 1.1.3
# Package for combining plots 

extend <- 
35*((0.25+(mahab[names]^.5)-min(mahab[names]^.5))/max(mahab[names]^.5))
# Scale mahalanobis distance for plotting 


data.mass <- data.frame(cbind(trait,mass,development[rownames(trait)])) 
colnames(data.mass)[c(1,3)]<-c('divergence','development')
data.mass[,1]<-as.numeric(data.mass[,1])
data.mass[,2]<-as.numeric(data.mass[,2])

linear.altricial.mass <-lm(mass~divergence,data=data.mass[which(data.mass$development=='Altricial'),])
ribbon.altricial.df.mass <-(seq(min(trait[names(which(development=='Altricial')),]),max(trait[names(which(development=='Altricial')),]),length.out=80))
ribbon.altricial.df.mass <- data.frame(ribbon.altricial.df.mass) ; colnames(ribbon.altricial.df.mass)<-'divergence' ; axis.altricial.mass <- data.frame(ribbon.altricial.df.mass[,1]); colnames(axis.altricial.mass)<-'divergence'
ribbon.altricial.df.mass[,2]<-  predict(linear.altricial.mass, newdata = axis.altricial.mass, se.fit = TRUE)$fit
ribbon.altricial.df.mass[,3]<- ribbon.altricial.df.mass[,2]- predict(linear.altricial.mass , newdata = axis.altricial.mass, se.fit = TRUE)$se.fit*1.96
ribbon.altricial.df.mass[,4]<- ribbon.altricial.df.mass[,2]+ predict(linear.altricial.mass , newdata = axis.altricial.mass, se.fit = TRUE)$se.fit*1.96
colnames(ribbon.altricial.df.mass)[2:4]<-c('y','ymin','ymax')
# summary(linear.altricial.mass)
# p < 0.001 divergence depending on mass in altricials for adjusted R^2 of 0.33 

linear.precocial.mass <-lm(mass~divergence,data=data.mass[which(data.mass$development=='Precocial'),])
ribbon.precocial.df.mass <-(seq(min(trait[names(which(development=='Precocial')),]),max(trait[names(which(development=='Precocial')),]),length.out=80))
ribbon.precocial.df.mass <- data.frame(ribbon.precocial.df.mass) ; colnames(ribbon.precocial.df.mass)<-'divergence' ; axis.precocial.mass <- data.frame(ribbon.precocial.df.mass[,1]); colnames(axis.precocial.mass)<-'divergence'
ribbon.precocial.df.mass[,2]<-  predict(linear.precocial.mass, newdata = axis.precocial.mass, se.fit = TRUE)$fit
ribbon.precocial.df.mass[,3]<- ribbon.precocial.df.mass[,2]- predict(linear.precocial.mass, newdata = axis.precocial.mass, se.fit = TRUE)$se.fit*1.96
ribbon.precocial.df.mass[,4]<- ribbon.precocial.df.mass[,2]+ predict(linear.precocial.mass, newdata = axis.precocial.mass, se.fit = TRUE)$se.fit*1.96
colnames(ribbon.precocial.df.mass)[2:4]<-c('y','ymin','ymax')
#ribbon.precocial.df[,2:4]<- 0.25+(ribbon.precocial.df[,2:4]- min(mahab[rownames(flight.pcoa$points)]^.5))
#ribbon.precocial.df[,2:4]<- 35*(ribbon.precocial.df[,2:4]/max(mahab[rownames(flight.pcoa$points)]^.5))
# We hypothesize that the rooted Mahalanobis distance will increase as a quadratic function of divergent wing and leg evolution; this 
# undertakes such a test in altricial birds and produces plotting metadata 
# summary(linear.precocial.mass)
# p ~ 0.56 divergence depending on mass in precocialss for adjusted R^2 of ~ 0


ggplot()+
	coord_flip()+
	labs( x='', y='')+
	geom_segment(data=phenogram[complete.cases(phenogram),],
  	aes(x = -y, y = -x, xend = -yend, yend =-xend), lwd = lwd.factor/10,col='yellow',lineend='round')+ 
	geom_point(aes(y=-trait[names,],x= (35*((0.25+(mahab[names]^.5)-min(mahab[names]^.5))/max(mahab[names]^.5))+1 ),
	fill= development[names],col= development[names]  ), shape=21, size=4*lwd.line,stroke=lwd.line*2 )+
	scale_fill_manual(values=c('#E69F00','#56B4E9'))+
	scale_color_manual(values=c('yellow','blue'))+
	geom_point(aes(y=-trait[names,][which(development[names]=='Precocial')],x= (35*((0.25+(mahab[names][which(development[names]=='Precocial')]^.5)-min(mahab[names]^.5))/max(mahab[names]^.5))+1 )),
	fill= '#56B4E9',col= 'blue'  , shape=21, size=4*lwd.line,stroke=lwd.line*2 )+

 	#geom_text(aes(y= jitter(-trait[names,], amount=.0) ,x= (35*((0.25+(mahab[names]^.5)-min(mahab[names]^.5))/max(mahab[names]^.5))+1 )),
	#label=gsub('_.*','',names))+
	# Should you wish to label the genus names of the species with described flight-styles. 

	geom_ribbon(data=ribbon.altricial.df,aes(y=trait,xmin=ymin+1,xmax=ymax+1),color=NA,fill='yellow',alpha=.4)+
	geom_path(data=ribbon.altricial.df,aes(y=trait,x=y+1),color='#E69F00',lwd=lwd.line,linetype='solid')+
	geom_ribbon(data=ribbon.precocial.df,aes(y=trait,xmin=ymin+1,xmax=ymax+1),color=NA,fill='#56B4E9',alpha=.2)+
	geom_path(data=ribbon.precocial.df,aes(y=trait,x=y+1),color='blue',lwd=lwd.line,linetype='dashed')+
	geom_segment(data=phenogram[complete.cases(phenogram),][which(phenogram[complete.cases(phenogram),]$col=='Altricial'),],
  	aes(x = -y, y = -x, xend = -yend, yend =-xend), lwd =lwd.factor/40 ,col='#E69F00',lineend='round')+
	geom_segment(data=phenogram[complete.cases(phenogram),][which(phenogram[complete.cases(phenogram),]$col=='Precocial'),],
  	aes(x = -y, y = -x, xend = -yend, yend =-xend), lwd = lwd.factor/10,col='blue',lineend='round')+ 
	geom_segment(data=phenogram[complete.cases(phenogram),][which(phenogram[complete.cases(phenogram),]$col=='Precocial'),],
  	aes(x = -y, y = -x, xend = -yend, yend =-xend), lwd = lwd.factor/40,col= '#56B4E9',lineend='round')+ 

	geom_point(aes(y=-trait[names(which(development[shortest$tip]=='Altricial')),],x=10*mass[names(which(development[shortest$tip]=='Altricial'))]-6+30), shape=22, size=3*lwd.line,stroke=1*lwd.line , show.legend = FALSE, fill='#E69F00',
	col='yellow')+

	geom_point(aes(y=-trait[names(which(development[shortest$tip]=='Precocial')),],x=10*mass[names(which(development[shortest$tip]=='Precocial'))]-6+30), shape=22, size=3*lwd.line,stroke=1*lwd.line , show.legend = FALSE, fill='#56B4E9',
	col='blue')+

	geom_ribbon(data=ribbon.altricial.df.mass,aes(y=-divergence,xmin=ymin*10-6+30,xmax=ymax*10-6+30),color=NA,fill='yellow',alpha=.4)+
	geom_path(data=ribbon.altricial.df.mass,aes(y=-divergence,x=y*10-6+30),color='#E69F00',lwd=lwd.line,linetype='solid')+
	geom_ribbon(data=ribbon.precocial.df.mass,aes(y=-divergence,xmin=ymin*10-6+30,xmax=ymax*10-6+30),color=NA,fill='#56B4E9',alpha=.2)+
	geom_path(data=ribbon.precocial.df.mass,aes(y=-divergence,x=y*10-6+30),color='blue',lwd=lwd.line,linetype='dashed')+
	geom_path(aes(x=c(0,35)+30,y=c(-0.32,-0.32)),lwd=lwd.line)+
	geom_path(aes(x=c(4,4)+30,y=c(-0.32,-0.34)),lwd=lwd.line)+
	geom_path(aes(x=c(14,14)+30,y=c(-0.32,-0.34)),lwd=lwd.line)+
	geom_path(aes(x=c(24,24)+30,y=c(-0.32,-0.34)),lwd=lwd.line)+
	geom_path(aes(x=c(34,34)+30,y=c(-0.32,-0.34)),lwd=lwd.line)+
	geom_text(aes(x=c(4,14,24,34)+30,y=rep(-0.36,4)),label=c('1','2','3','4'), size=4)+
	geom_text(aes(x=c(38)+30,y=-0.3),label=expression(paste(italic('mass'),' (log'[10],' g)')), size=4)+


	draw_plot(bars, y = .4, x = -105, height = .45, width = 65)+ 

	geom_path(aes(x=c(0,max(extend)),y=c(0.5,0.5)),lwd=lwd.line)+
	geom_path(aes(x=c(max(extend)*0.25,max(extend)*0.25),y=c(0.5,0.52)),lwd=lwd.line)+
	geom_path(aes(x=c(max(extend)*0.5,max(extend)*0.5),y=c(0.5,0.52)),lwd=lwd.line)+
	geom_path(aes(x=c(max(extend)*0.75,max(extend)*0.75),y=c(0.5,0.52)),lwd=lwd.line)+
	geom_path(aes(x=c(max(extend),max(extend)),y=c(0.5,0.52)),lwd=lwd.line)+
	geom_text(aes(x=c(max(extend)*.25,max(extend)*.5,max(extend)*.75,max(extend)),y=rep(0.55,4)),label=c('0.25','0.5','0.75','1'),size=4)+
	geom_text(aes(x=c(31),y=0.5),label=expression(sqrt(italic(M[flight]))),size=4 )+
	guides(fill = guide_legend(override.aes = list(size=3)))+
	xlab('Mya')+
	ylab(expression(paste(hat(Delta)," (wing-leg divergence)")))+
	geom_text( aes(x=c(73.5,24,-20),y=c(-0.275,-0.275,-0.275)),label=c('(a)','(b)','(c)'), size = 11/.pt, color = "black", fontface = "bold.italic")+
	lims(x=c(-104,74))+
	theme(
	legend.position=c(0.12,0.20),
	legend.title=element_blank(),
	legend.key = element_rect(fill = alpha("white", 0.0)),
         panel.background = element_rect(fill='white'),
         plot.background = element_rect(fill='white', color=NA),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
		legend.background = element_rect(fill='transparent'),
		legend.key.width=unit( 0.4,'cm'),
		legend.key.height=unit( 0.5,'cm'),
		legend.text=element_text(size=11,color='black' ),
		axis.line.x=element_line(color='black',linewidth=lwd.line),
		axis.line.y=element_line(color='black',linewidth=lwd.line),
		axis.ticks.x=element_line(color='black',linewidth=lwd.line),
		axis.ticks.y=element_line(color='black',linewidth=lwd.line),
		axis.title.x=element_text(size=11, color='black'),
		axis.title.y=element_text(size=11, color='black'),
		axis.text.x=element_text(size=11,color='black'),
		axis.text.y=element_text(size=11,color='black'),

)


ggsave(filename='Figure_2_07_31_2025.pdf',height=18,width=18,unit='cm', dpi=300,device="pdf",family='Times') 

# This visualization makes it apparent that Precocial birds diverge less from a concordance axis of wing and leg evolution
# Altricial birds explore divergence from this axis, including forms with relatively very small or very large wings.
# Unusual divergences between wing and leg evolution in altricial birds are signifcantly associated with unusual 
# combinations of flight style characteristics.
# The birds with the most unusually large wings are Frigate birds, tropicbirds, Pelecans and Gannets. They do interesting things like
# soaring over the ocean and pursuing other targets in the air (e.g. during piracy). 
# The birds with the most unusually short wings relative to their legs are Hummingbirds. They practice 
# interesting combinations of flight style such as Hovering, Pursuing and Sallying for insects 
# Birds which sally and flap-bound in cluttered environments are not very unusual. 
# Fairy wrens, owls, parrots and passerines tend to have fairly conservative limb proportions and flight styles. 
# Very unusual birds which have unique combinations of flight styles in spite of having conventional limb proportions
# are limited, such as Puffinus Griseus, which is unusual because it will dive tens of meters into the sea and propel
# itself with its wings. 
# Some birds, such as Pelagodroma_marina, probably have more unusual flight styles than what is represented here.
# This bird has very unusually long legs compares to its wings, and it dances across the sea surface in an action
# known as 'pattering', to feel for tiny crustaceans with its feet. This is not a scored flight character, but clearly
# shows this bird has an unusual relationship between aberrant limb evolution and unique foraging ecology 

# Notice that the Precocial birds in general do not obtain the most conventional flight style states;
# this is because the 'mean' bird is some passerine-like bird which lives in clutter and flap-bounds, 
# which tend to be less accessible to precocial birds, that are often tied to water environments or which 
# practice burst flight. 


# Script concludes. 
