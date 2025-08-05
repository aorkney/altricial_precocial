# The purpose of this script is to explore the strength of evolutionary covariance between 
# the relative proportions of wing and leg long bone lengths across birds. 
# We want to know whether these patterns differ substantially between bird lineages with high parental investment (altricial)
# and lineages with low parental investment (precocial). 
# We will present these results as a heatmap illustrating the evolutionary correlations between 
# different components of the avian limb skeleton in altricial and precocial birds. 
# 
# Patterns of evolutionary covariance, or 'modules' are of interest to biologists, 
# because they may represent the manifestation of restrictions on the correlated growth, development or selection 
# across groups of traits and therefore may either restrain or permit different sorts of adaptation. 
#
# We hypothesise that the evolutionary covariance between wing and leg bones will be weaker in altricial birds


library(ggplot2) # 3.5.1
# Package for producing plots

# Specify some plot parameters 
 Subplot.title.font.size <- 10 
 Subplot.axis.font.size <- 9 
 Legend.axis.font.size <- 7.5 
 Legend.title.font.size <- 9 
 Legend.tick.size <- 1/4
 Subplot.module.linewidth <- 1/4
 Subplot.tick.size <-1/4
 Subplot.border.linewidth <- 1/4
 legend.key.width <- 1/2
 legend.key.height <- 0.5

library( ape ) # 5.8-1 
library(phytools) # 2.4-4
# Load required analytical R packages.


setwd('')
# Set work directory


library(clootl) # 0.1.1
# The McTavish et al. 2025 avian super tree will be used to represent
# evolutionary relationships between birds:
# https://doi.org/10.1073/pnas.2409658122

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


development <- rep(NA,dim(combined)[1])
development[grep('tricial_1|tricial_2',combined$Ontogeny)]<-'Altricial'
development[grep('recocial_1|recocial_2|recocial_3|uperrecocial',combined$Ontogeny)]<-'Precocial'
names(development)<-combined$Species
development<-development[short.tree$tip.label]
# Define a vector delineating birds which are strongly precocial or altricial in their developmental style
# Sourced from Starck 1993

pairs <- combn(colnames(r),2)
# Define pairs of bones to compute the evolutionary covariance between

precocial.battenberg<- matrix(NA,dim(r)[2],dim(r)[2])
rownames(precocial.battenberg)<-colnames(precocial.battenberg)<-colnames(r)
altricial.battenberg<- precocial.battenberg
diff.p<- altricial.battenberg
# Define matrices to store the results of our analyses


# Compute pearson correlations of residual bone lengths in an evolutionary context
# We will also perform a significance test to determine if they are distinct between
# Altricial and Precocial cohorts of birds:

# Compute pearson correlations of transformed data pairs
cor.score <- function(r, phy, group1,group2){
	C <- ape::vcv.phylo(phy)
	r<- r[rownames(C),c(group1,group2)]
	# Order data
	eigen <- eigen(C)
	# Eigendecomposition
	invC <- solve(C)
	# Find inverse of C
	p.tran <- solve(eigen$vectors %*% diag(sqrt(eigen$values)) %*% t(eigen$vectors))
	p.tran.r <- p.tran %*% r
	# Transform trait data into basis aligned with shared ancestry

	XScores <- p.tran.r[,1] 
	YScores <- p.tran.r[,2] 
	# Projected traits

	return(
	list(dim(C)[1],cor(XScores,YScores),
	sd(XScores)*sd(YScores)	) )
	# Return the Pearson's correlation of the transformed data and other test statistics for downstream use

}



	C <- ape::vcv.phylo(keep.tip(short.tree,names(which(development=='Altricial'))))
	j<-1
	Vr <- 0
	for(j in 1:dim(C)[1]){
		Vr <- Vr + (1-C[j,-j]%*%(solve(C[-j,-j])/(C[j,j] ))%*%C[-j,j])
	}
	Neff.A <- 1+ ((dim(C)[1]-1)/dim(C)[1]) * Vr
	# https://arxiv.org/pdf/1507.07113
	# Impute the effective sample size, accommodating shared ancestry.
	# Note that I am making a simplifying assumption that shared ancestry represents 
	# variance-covariance in traits adequately as a Brownian process



	C <- ape::vcv.phylo(keep.tip(short.tree,names(which(development=='Precocial'))))
	j<-1
	Vr <- 0
	for(j in 1:dim(C)[1]){
		Vr <- Vr + (1-C[j,-j]%*%(solve(C[-j,-j])/(C[j,j] ))%*%C[-j,j])
	}
	Neff.P <- 1+ ((dim(C)[1]-1)/dim(C)[1]) * Vr
	# https://arxiv.org/pdf/1507.07113
	# Impute the effective sample size, accommodating shared ancestry.
	# Note that I am making a simplifying assumption that shared ancestry represents 
	# variance-covariance in traits adequately as a Brownian process


# The following set of loops applies our function across all pairwise combinations of bones
# for precocial and altricial birds, and explicitly uses Fisher transformed correlation values
# and effective sample size estimates to test if the correlations are significantly different
# between altricial and precocial cohorts. 
# Specifically, we test the one-tailed hypothesis that the correlations are weaker in altricials. 
# Be prepared for this loop to take perhaps 30 minutes: 

for(i in 1:dim(pairs)[2]){
	bone1 <- pairs[1,i]
	bone2 <- pairs[2,i]
	
	c.precocial <- cor.score(r,phy=keep.tip(short.tree,names(which(development=='Precocial')) ),group1=bone1,group2=bone2)
	c.altricial <- cor.score(r,keep.tip(short.tree,names(which(development=='Altricial')) ),group1=bone1,group2=bone2)
	z <- (atanh(c.precocial[[2]])- atanh(c.altricial[[2]]) )/
	sqrt( (1/(Neff.P)) + (1/(Neff.A)) )
	pnorm(z,lower.tail=F)


	row <- which(rownames(precocial.battenberg)==bone1)
	column <- which(colnames(precocial.battenberg)==bone2)
	precocial.battenberg[row,column] <- precocial.battenberg[column,row] <- c.precocial[[2]] 
	altricial.battenberg[row,column] <- altricial.battenberg[column,row] <- c.altricial[[2]] 
	diff.p[row,column] <- diff.p[column,row] <- pnorm(z,lower.tail=F)
}



library(reshape2) # 1.4.4
# This package contains functions for formatting data 
pdf <- melt(precocial.battenberg)
colnames(pdf)<-c('bone1','bone2','Z')
# Prepare precocial statistical results

adf <- melt(altricial.battenberg)
colnames(adf)<-c('bone1','bone2','Z')
# Prepare altricial statistical results

ddf <- melt(diff.p)


limits.one <- c(0,1)
# Limits of valid correlation coefficients

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# This is a colour blind friendly palette. 

bone_colours<-rev(c(rep(cbbPalette[4],3),rep(cbbPalette[8],3)))
# Colours for different modules.

# Heatmap for precocial birds
pre <- 
ggplot(data=pdf)+
geom_tile(aes(x=factor(bone1, levels= colnames(r)), y=factor(bone2, levels= colnames(r)), fill=Z))+
geom_point(data=pdf[which(ddf$val<0.05),], aes(x=factor(bone1, levels= colnames(r)), y=factor(bone2, levels= colnames(r))),shape="asterisk",size=4,stroke=2,color='white')+
scale_y_discrete( limits= rev,labels=substr(rev(colnames(r)),1,2) )+
scale_x_discrete(labels=substr(colnames(r),1,2), guide = guide_axis(angle = 90))+
scale_fill_gradient2(high='black',low='white',midpoint=0.0,na.value='white', limits=limits.one )+
theme(axis.text.x=element_text(size=Subplot.axis.font.size, angle=90, vjust=0.3,colour=rev(bone_colours),face = "bold"),
                     axis.text.y=element_text(colour=(bone_colours),face = "bold"),
				#axis.text.y=element_blank(),
				axis.title.x=element_blank(),
 				axis.ticks=element_line(size=1/2),
				axis.title.y=element_blank(),
				legend.key.width=unit(legend.key.width*2,'cm'),
				legend.key.height=unit(legend.key.height,'cm'),  
				legend.title=element_text(size=Legend.title.font.size*2),
				legend.text=element_text(size=Legend.axis.font.size),
				plot.margin=unit(c(0.1,0.1,0.1,0.1),'cm'),
                     plot.title=element_text(size=Subplot.title.font.size,hjust=.5), legend.position="bottom",panel.background = element_blank())+
	labs(fill = expression(italic(rho)) )+
geom_segment(data=as.data.frame(cbind(x=c(0.5,0.5,0.5,0.5,3.5,6.5),xend=c(6.5,6.5,6.5,0.5,3.5,6.5),y=c(0.5,3.5,6.5,0.5,0.5,0.5),yend=c(0.5,3.5,6.5,6.5,6.5,6.5) )), 
aes(x=x,xend=xend,y=y,yend=yend), size=1, colour='black')+
annotate('rect',colour='black',fill=NA,xmin=3.5,xmax=6.5,ymin=0.5,ymax=3.5,linewidth=Subplot.module.linewidth+1/4)+
annotate('rect',colour=cbbPalette[8],fill=NA,xmin=3.5,xmax=6.5,ymin=0.5,ymax=3.5,linewidth=Subplot.module.linewidth)+
annotate('rect',colour='black',fill=NA,xmin=0.5,xmax=3.5,ymin=3.5,ymax=6.5,linewidth=Subplot.module.linewidth)+
annotate('rect',colour=cbbPalette[4],fill=NA,xmin=0.5,xmax=3.5,ymin=3.5,ymax=6.5,linewidth=Subplot.module.linewidth+1/4)+
ggtitle(expression(paste(italic('(a)'),' Precocial' ,italic(' n'),' = 222')))+
 guides(fill = guide_colourbar( ticks.linewidth=Legend.tick.size, frame.colour = 'white',
  frame.linewidth = 0/.pt))+ coord_fixed()


battenberg_legend <- cowplot::get_legend(pre)
# Extract the legend for future plotting. 
battenberg_legend <-cowplot::get_plot_component(pre,pattern='guide-box',return_all = TRUE)[[3]]

# Heatmap for precocial birds without legend
pre <- 
ggplot(data=pdf)+
geom_tile(aes(x=factor(bone1, levels= colnames(r)), y=factor(bone2, levels= colnames(r)), fill=Z))+
geom_point(data=pdf[which(ddf$val<0.05),], aes(x=factor(bone1, levels= colnames(r)), y=factor(bone2, levels= colnames(r))),shape="asterisk",size=1,stroke=1/2,color='white')+
scale_y_discrete( limits= rev,labels=substr(rev(colnames(r)),1,2) )+
scale_x_discrete(labels=substr(colnames(r),1,2), guide = guide_axis(angle = 90))+
scale_fill_gradient2(high='black',low='white',midpoint=0.0,na.value='white', limits=limits.one )+
theme(plot.background = element_rect(fill='transparent', color=NA), axis.text.x=element_text(size=Subplot.axis.font.size, angle=90, vjust=0.3,colour=rev(bone_colours),face = "bold"),
                     axis.text.y=element_text(colour=(bone_colours),face = "bold"),
				#axis.text.y=element_blank(),
				axis.title.x=element_blank(),
 				axis.ticks=element_line(size=1/2),
				axis.title.y=element_blank(),
				legend.key.width=unit(legend.key.width,'cm'),
				legend.key.height=unit(legend.key.height,'cm'),  
				legend.title=element_text(size=Legend.title.font.size),
				legend.text=element_text(size=Legend.axis.font.size),
				plot.margin=unit(c(0.1,0.1,0.1,0.1),'cm'),
                     plot.title=element_text(colour= "#56B4E9",size=Subplot.title.font.size,hjust=.5), legend.position="none",panel.background = element_blank())+
	labs(fill = expression(italic(rho)) )+
geom_segment(data=as.data.frame(cbind(x=c(0.5,0.5,0.5,0.5,3.5,6.5),xend=c(6.5,6.5,6.5,0.5,3.5,6.5),y=c(0.5,3.5,6.5,0.5,0.5,0.5),yend=c(0.5,3.5,6.5,6.5,6.5,6.5) )), 
aes(x=x,xend=xend,y=y,yend=yend), size=1/4, colour='black')+
annotate('rect',colour='black',fill=NA,xmin=3.5,xmax=6.5,ymin=0.5,ymax=3.5,linewidth=Subplot.module.linewidth+1/4)+
annotate('rect',colour=cbbPalette[8],fill=NA,xmin=3.5,xmax=6.5,ymin=0.5,ymax=3.5,linewidth=Subplot.module.linewidth)+
annotate('rect',colour='black',fill=NA,xmin=0.5,xmax=3.5,ymin=3.5,ymax=6.5,linewidth=Subplot.module.linewidth)+
annotate('rect',colour=cbbPalette[4],fill=NA,xmin=0.5,xmax=3.5,ymin=3.5,ymax=6.5,linewidth=Subplot.module.linewidth+1/4)+
ggtitle(expression(paste(italic('(a)'),' Precocial' ,italic(' n'),' = 222')))+
 guides(fill = guide_colourbar( ticks.linewidth=Legend.tick.size, frame.colour = 'white',
  frame.linewidth = 0/.pt))+ coord_fixed()

# Heatmap for altricial birds
alt <- 
ggplot(data=adf)+
geom_tile(aes(x=factor(bone1, levels= colnames(r)), y=factor(bone2, levels= colnames(r)), fill=Z))+
geom_point(data=adf[which(ddf$val<0.05),], aes(x=factor(bone1, levels= colnames(r)), y=factor(bone2, levels= colnames(r))),shape="asterisk",size=1,stroke=1/2,
color=c('white','black','black','white','black','black','black','black','black','black','black','black','black','black','black','white','black','black','black','white'))+
scale_y_discrete( limits= rev,labels=substr(rev(colnames(r)),1,2) )+
scale_x_discrete(labels=substr(colnames(r),1,2), guide = guide_axis(angle = 90))+
scale_fill_gradient2(high='black',low='white',midpoint=0.0,na.value='white', limits=limits.one )+
theme(plot.background = element_rect(fill='transparent', color=NA), axis.text.x=element_text(size=Subplot.axis.font.size, angle=90, vjust=0.3,colour=rev(bone_colours),face = "bold"),
                     axis.text.y=element_text(colour=(bone_colours),face = "bold"),
				#axis.text.y=element_blank(),
				axis.title.x=element_blank(),
 				axis.ticks=element_line(size=1/2),
				axis.title.y=element_blank(),
				legend.key.width=unit(legend.key.width,'cm'),
				legend.key.height=unit(legend.key.height,'cm'),  
				legend.title=element_text(size=Legend.title.font.size),
				legend.text=element_text(size=Legend.axis.font.size),
				plot.margin=unit(c(0.1,0.1,0.1,0.1),'cm'),
                     plot.title=element_text(colour= '#E69F00',size=Subplot.title.font.size,hjust=.5), legend.position="none",panel.background = element_blank())+
	labs(fill = expression(italic(rho)) )+
geom_segment(data=as.data.frame(cbind(x=c(0.5,0.5,0.5,0.5,3.5,6.5),xend=c(6.5,6.5,6.5,0.5,3.5,6.5),y=c(0.5,3.5,6.5,0.5,0.5,0.5),yend=c(0.5,3.5,6.5,6.5,6.5,6.5) )), 
aes(x=x,xend=xend,y=y,yend=yend), size=1/4, colour='black')+
annotate('rect',colour='black',fill=NA,xmin=3.5,xmax=6.5,ymin=0.5,ymax=3.5,linewidth=Subplot.module.linewidth+1/4)+
annotate('rect',colour=cbbPalette[8],fill=NA,xmin=3.5,xmax=6.5,ymin=0.5,ymax=3.5,linewidth=Subplot.module.linewidth)+
annotate('rect',colour='black',fill=NA,xmin=0.5,xmax=3.5,ymin=3.5,ymax=6.5,linewidth=Subplot.module.linewidth)+
annotate('rect',colour=cbbPalette[4],fill=NA,xmin=0.5,xmax=3.5,ymin=3.5,ymax=6.5,linewidth=Subplot.module.linewidth+1/4)+
ggtitle(expression(paste(italic('(b)'),' Altricial' ,italic(' n'),' = 386')))+
 guides(fill = guide_colourbar( ticks.linewidth=Legend.tick.size, frame.colour = 'white',
  frame.linewidth = 0/.pt))+ coord_fixed()

# Note that asterisks indicate that correlations are significantly lower in altricial birds 

# Now it is time to load illustrated visual elements (not included for the public version)
#library(jpeg)
# Package for processing jpeg image files

#precocial.img <- readJPEG("precocial_bird.jpeg")
#altricial.img <- readJPEG("altricial_bird.jpeg")
# Illustrations loaded

library(cowplot)
# Useful package for combining images and subplots

blank<-ggplot()+ geom_blank()+
  theme(
	axis.text.x=element_text(size=12,colour='black'),
  	axis.text.y=element_text(size=12,colour='black'),
	axis.title.x.bottom=element_text(size=15,colour="black"),
	axis.title.y.left=element_text(size=15,colour="black"),
      axis.ticks=element_line(size=2),
      panel.background=element_blank(),
      panel.border=element_blank(),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      plot.background=element_blank())
# Done


  ggdraw() +
  draw_plot(blank) +
 # draw_image(altricial.img,x=0.37,y=0.10,width=0.8,height=.85)+
 # draw_image(precocial.img,x=-0.245,y=0.055,width=0.8,height=.85)+
  draw_plot(pre, x = 0.3625, y = 0.0925, width = .2, height = 1)+
  draw_plot(alt, x = 0.3625, y = -0.28, width = .2, height = 1)+ #-0.175
	draw_plot(battenberg_legend,x=0.375,y=0.7,width=0.2,height=0.35)
 # Total compiled plot


ggsave(filename='Figure_1_07_31_2025.pdf',width=18,height=11,unit='cm',dpi=600,device='pdf',family='Times') 
# Save the plot

# Script concludes 



