# The purpose of this script is to produce a ternary diagram representing
# internal wing proportions across a representative sample of birds.
# We will visualise how the magnitude of evolutionary divergence between
# wings and legs is distributed across this space, 
# how a discretised gradient of precocial (low parental care) to altricial
# (high parental care) is distributed, and subset
# the data to several large avian clades to check that perceived patterns
# are general. 

setwd('')
# Set work directory to the location that you downloaded the data 

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


reorient <- solve(svd(cbind(XScores[,1],YScores[,1]))$v)
# Rotate the distribution of accorded wing-leg evolution so that the second axis described departure from concordance
# reorient %*% t(cbind(pls$u[,1],pls$v[,1])) # Inspect reoriented coefficients if desired 
#           [,1]       [,2]        [,3]
#[1,]  0.7752074  0.7883293  0.84396839
#[2,] -0.1674712 -0.1724931 -0.08664417
# The reoriented coefficients demonstrate that the concordance axis is dominated by correlated
# evolution between the wing and more proximal components of the leg, and in particular that the 
# zeugopodial component of the wing (radius) is strongly related to the leg. 

plot(t(reorient %*% t(cbind(XScores[,1],YScores[,1]))),lwd=4)

divergence <- t(reorient %*% t(cbind(XScores[,1],YScores[,1])))
rownames(divergence) <- short.tree$tip
trait <- as.matrix(divergence[,2])
# Define our axis of divergent wing and leg evolution


development <- combined$Ontogeny
names(development)<-combined$Species
development<-development[short.tree$tip]

lengths <- combined[match(short.tree$tip.label,combined$Species),c(6,8,9,10,11,12)]
rownames(lengths) <- short.tree$tip.label

df<-data.frame(lengths)
colnames(df)<-c('humerus','radius','carpometacarpus','femur','tibiotarsus','tarsometatarsus')
df<-apply(df,2,as.numeric)
df[,1:3]<-df[,1:3]/rowSums(df[,1:3])
df[,4:6]<-df[,4:6]/rowSums(df[,4:6])

df<-cbind(df,development[rownames(lengths)] )
colnames(df)[7]<-'mode'

df<-as.data.frame(df)
# Bind the data into a dataframe and add developmental strategy categorisations. 
df$mode<- factor(df$mode, levels=c('Superprecocial', 'Precocial_1', 'Precocial_2', 'Precocial_3', 'Semiprecocial', 'Semialtricial', 'Altricial_1', 'Altricial_2') )

df<- cbind(df, combined$Mag[match(rownames(lengths),combined$Species)])
colnames(df)[8] <-'clade'

# We need to take action here to complete the clade labels 
nas <- which(is.na(df$clade)==T)
# incomplete labels 

df$clade[nas][which(nas<=289)]<-'Telluraves'
df$clade[nas][which(nas>289 & nas <310)]<-'Strisores'
df$clade[nas][which(nas>310 & nas <433)]<-'Phaethoquornithes'
df$clade[nas][which(nas>433 & nas <558)]<-'Cursorimorphae'
df$clade[nas][which(nas>559 & nas <631)]<-'Columbaves'
df$clade[nas][which(nas>642 & nas <764)]<-'Galloanseres'
df$clade[nas][which(nas>763)]<-'Palaeognathae'
# completed 

df$humerus <- as.numeric(df$humerus)
df$radius <- as.numeric(df$radius)
df$carpometacarpus <- as.numeric(df$carpometacarpus)
df$femur <- as.numeric(df$femur)
df$tibiotarsus <- as.numeric(df$tibiotarsus)
df$tarsometatarsus <- as.numeric(df$tarsometatarsus)

# Function to convert degrees to radians 
deg2rad <- function(deg) {(deg * pi) / (180)}

# Function for trigonometry for producing ternary diagrams
tern2cart <- function(m1,m2){
	x <- -(m1 + (m2)*cos(deg2rad(60)))
	y <- m2*sin(deg2rad(60))
	return(cbind(x,y))
}

dfh <-data.frame(tern2cart(df$humerus*100,df$radius*100))
dfh <- cbind(dfh, df[,7:8])
# Do the trigonometry 

# Plot aesthetics 
outline <- data.frame(tern2cart(c(20,20, 40,60,20),c(30,50, 50,30,30)))
lines <- rbind(cbind(tern2cart(rep(20,9),seq(30,50,length.out=9)), tern2cart(seq(60,40,length.out=9),seq(30,50,length.out=9)) ),
cbind(tern2cart(seq(20,60,length.out=17),rep(30,17)), tern2cart(seq(20,60,length.out=17),c(rep(50,8),seq(50,30,length.out=9) ) ) ),
cbind(tern2cart(seq(60,20,length.out=17),rep(30,17)) , tern2cart(c(seq(40,20,length.out=9),rep(20,8)),c(rep(50,8),seq(50,30,length.out=9))) ) )
colnames(lines)<-c('X1','X2','X3','X4')
ticks <- rbind( cbind(tern2cart(rep(20,9),seq(30,50,length.out=9)),tern2cart(rep(19,9),seq(30,50,length.out=9))),
cbind(tern2cart(seq(20,60,length.out=17),rep(30,17)),tern2cart(seq(21,61,length.out=17),rep(30-1,17)) ),
cbind(tern2cart(c(seq(20,40,length.out=9),seq(40,60,length.out=9)) ,c(rep(50,9),seq(50,30,length.out=9))) ,tern2cart(c(seq(20,40,length.out=9),seq(40,60,length.out=9)) ,c(rep(51,9),seq(51,31,length.out=9))) )  )
colnames(ticks)<-c('X1','X2','X3','X4')

n<-100
zx <- log(df[,1:3]) # log transform the proportions
z <- as.matrix(zx-rowMeans(zx)) # center by subtracting means
m <- colMeans(z) # mean vector
a <- eigen(cov(z))$vectors[,1] + m # first unit vector, centered on mean
sc <- z %*% a # Are these the loadings? It is the projection of the data on the first unit vector
 lam <- seq( min(sc) - .1, max(sc) - .55, length = n )
    x1 <- cbind( a[1] * lam, a[2] * lam, a[3] * lam) + cbind( m[1] * (1 - lam),
    m[2] * (1 - lam), m[3] * (1 - lam) )
    expx1 <- exp(x1)
    wa1 <- expx1 /rowSums( expx1 )  ## first principal component in S^2
dfline <-data.frame(tern2cart(wa1[,1]*100,wa1[,2]*100))

zx <- log(df[,1:3]) # log transform the proportions
z <- as.matrix(zx-rowMeans(zx)) # center by subtracting means
a <- eigen(cov(z))$vectors[,2] + m # first unit vector, centered on mean
sc <- z %*% a # Are these the loadings? It is the projection of the data on the first unit vector
 lam <- seq( min(sc) + .02, max(sc) - .05, length = n )
    x1 <- cbind( a[1] * lam, a[2] * lam, a[3] * lam) + cbind( m[1] * (1 - lam),
    m[2] * (1 - lam), m[3] * (1 - lam) )
    expx1 <- exp(x1)
    wa1 <- expx1 /rowSums( expx1 )  ## first principal component in S^2
dfline2 <-data.frame(tern2cart(wa1[,1]*100,wa1[,2]*100))

library(ggplot2) # 3.5.1


# Plot parameters 
fontsize =(11/.pt)/2 # , family = "Times New Roman"
line.width = 1/4
outline.width = 3/4
fineline.width = 1/8
point.size = 1
fill.size =1/2




dfh$divergence <- (abs(divergence[,2]-mean(divergence[,2]))^.5)


# Ternary diagram visualising dispersion of the relative magnitude of evolutionary divergence between the wing and leg. 

subplota <- 
ggplot( ) +
scale_x_continuous(limits=range(c(ticks[,1],ticks[,3]))+c(-5,5),expand=c(0,0))+
geom_segment(data=lines,aes(x=X1,xend=X3,y=X2,yend=X4),linewidth=fineline.width,col='grey',alpha=1)+
geom_path(data=outline,aes(x=x,y=y),linewidth=outline.width,col='black')+
labs(fill = expression(italic( bold("|")~hat(Delta)-bar(Delta)~bold("|") ^~frac(1,2) ))    )+
geom_segment(data=ticks,aes(x=X1,xend=X3,y=X2,yend=X4),linewidth=outline.width,alpha=1, col=  c(viridis::magma(13)[1:9],rev(viridis::mako(25)[5:21]),viridis::viridis(22)[1:18]) )+
geom_path(data=dfline,aes(x=x,y=y),size=outline.width,color='black',arrow = arrow(angle = 15, ends = "both", type = "closed",length=unit(0.2,'cm')))+
geom_path(data=dfline2,aes(x=x,y=y),size=outline.width,color='black',arrow = arrow(angle = 15, ends = "both", type = "closed",length=unit(0.2,'cm')))+
geom_point(data=dfh,aes(x=x,y=y),size=point.size,color='black')+
geom_point(data=dfh[rev(order(dfh$divergence)),],aes(x=x,y=y, fill= divergence ),size=fill.size,color='NA',shape=21)+
scale_fill_gradientn(colors=rainbow(100)[45:90])+
geom_text( aes(x= -(c(18.5,18.5) + c(30,50)*cos(deg2rad(60))),y= c(30,50)*sin(deg2rad(60)) ) ,label=c(30,50),size=fontsize )+
geom_text( aes(x= -(c(55,25) + c(30,30)*cos(deg2rad(60))),y= c(28.5,28.5)*sin(deg2rad(60)) ) ,label=c(15,45),angle=240 ,size=fontsize)+
geom_text( aes(x= -(c(21,41,61) + c(50,50,30)*cos(deg2rad(60))),y= c(51,51,31)*sin(deg2rad(60)) ) ,label=c(20,40,60),angle=120,size=fontsize )+
geom_text( aes(x= -(c(18.5,17.5,60) + c(50,32,32)*cos(deg2rad(60))),y= c(52,34,33)*sin(deg2rad(60)) ) ,label=c('Radius','Carpometacarpus','Humerus'),angle=c(0,240,120), color=c(viridis::magma(13)[8],viridis::mako(25)[16],viridis::viridis(22)[16]) ,size=fontsize )+
geom_text(aes(x= -(32.5 + 35*cos(deg2rad(60))), y= 35*sin(deg2rad(60)) ) , label='PC1',angle=330,size=fontsize)+
geom_text(aes(x= -(32.5 + 46*cos(deg2rad(60))), y= 46*sin(deg2rad(60)) ) , label='PC2',angle=60,color='black',size=fontsize)+
#geom_text(data=dfh,aes(x=x,y=y),label=round(df[rownames(dfh),3],digits=2) )+
theme( plot.margin = unit(c(-10,0,-10,0), 'lines'),   
 panel.background = element_rect(fill='transparent'),
      plot.background = element_rect(fill='transparent', color=NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
	axis.title.x=element_blank(),
	axis.title.y=element_blank(),
	axis.text.x=element_blank(),
	axis.text.y=element_blank(),
	axis.ticks.y=element_blank(),
	axis.ticks.x=element_blank(),
	legend.position=c(0.1,0.65),
	legend.text=element_text(size=fontsize*.pt),
	legend.title=element_text(size=fontsize*.pt),
	legend.direction="vertical",
	legend.background = element_rect(fill = NA) )+
 coord_fixed()


df$mode<- factor(df$mode, levels=c('Superprecocial', 'Precocial_1', 'Precocial_2', 'Precocial_3', 'Semiprecocial', 'Semialtricial', 'Altricial_1', 'Altricial_2') )
# Stack the order of the discretised gradient of developmental strategies

# Ternary diagram visualising the developmental strategies across the ternary space: 
subplotb <- 
ggplot( ) +
geom_segment(data=lines,aes(x=X1,xend=X3,y=X2,yend=X4),linewidth=fineline.width,col='grey',alpha=1)+
geom_path(data=outline,aes(x=x,y=y),linewidth=outline.width,col='black')+
labs(fill = ''    )+
geom_segment(data=ticks,aes(x=X1,xend=X3,y=X2,yend=X4),linewidth=outline.width,alpha=1, col=  c(viridis::magma(13)[1:9],rev(viridis::mako(25)[5:21]),viridis::viridis(22)[1:18]) )+
geom_path(data=dfline,aes(x=x,y=y),size=outline.width,color='black',arrow = arrow(angle = 15, ends = "both", type = "closed",length=unit(0.2,'cm')))+
geom_path(data=dfline2,aes(x=x,y=y),size=outline.width,color='black',arrow = arrow(angle = 15, ends = "both", type = "closed",length=unit(0.2,'cm')))+geom_point(data=dfh,aes(x=x,y=y),size=point.size,color='black')+
geom_point(data=dfh[rev(order(dfh$mode)),],aes(x=x,y=y, fill=factor(mode) ),size=fill.size,color='NA',shape=21)+
scale_discrete_manual(breaks=unique(df$mode)[c(8,4,7,3,6,2,5,1)], aesthetics='fill',values=c(rev(rainbow(100)[c(1,10,15,65,60,58,56,50)]))[c(5,1,4,6,3,7,2,8)] )+
geom_point(data=dfh['Macrocephalon_maleo',],aes(x=x,y=y ),fill=rev(rainbow(100)[c(1,10,15,65,60,58,56,50)])[5], size=fill.size*6,color='white',shape=21)+
geom_text( aes(x= -(c(18.5,18.5) + c(30,50)*cos(deg2rad(60))),y= c(30,50)*sin(deg2rad(60)) ) ,label=c(30,50), size=fontsize )+
geom_text( aes(x= -(c(55,25) + c(30,30)*cos(deg2rad(60))),y= c(28.5,28.5)*sin(deg2rad(60)) ) ,label=c(15,45),angle=240,size=fontsize)+
geom_text( aes(x= -(c(21,41,61) + c(50,50,30)*cos(deg2rad(60))),y= c(51,51,31)*sin(deg2rad(60)) ) ,label=c(20,40,60),angle=120,size=fontsize )+
geom_text( aes(x= -(c(18.5,17.5,60) + c(50,32,32)*cos(deg2rad(60))),y= c(52,34,33)*sin(deg2rad(60)) ) ,label=c('Radius','Carpometacarpus','Humerus'),angle=c(0,240,120), color=c(viridis::magma(13)[8],viridis::mako(25)[16],viridis::viridis(22)[16]) ,size=fontsize )+
geom_text(aes(x= -(32.5 + 35*cos(deg2rad(60))), y= 35*sin(deg2rad(60)) ) , label='PC1',angle=330,size=fontsize)+
geom_text(aes(x= -(32.5 + 46*cos(deg2rad(60))), y= 46*sin(deg2rad(60)) ) , label='PC2',angle=60,color='black',size=fontsize)+
guides(fill = guide_legend(label.position = "bottom"))+
guides(fill = guide_legend(override.aes = list(size=4)))+
theme( plot.margin = unit(c(-10,0,-10,0), 'lines'),   
 panel.background = element_rect(fill='transparent'),
      plot.background = element_rect(fill='transparent', color=NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
	axis.title.x=element_blank(),
	axis.title.y=element_blank(),
	axis.text.x=element_blank(),
	axis.text.y=element_blank(),
	axis.ticks.y=element_blank(),
	axis.ticks.x=element_blank(),
	legend.position=c(0.5,1.1),
	legend.text=element_text(size=fontsize*.pt,angle=0),
	legend.title=element_text(size=fontsize*.pt),
	legend.direction="horizontal",
	legend.spacing.x = unit(0, 'cm'),
	legend.background = element_rect(fill = NA)  )+
 coord_fixed()


#library(png) # 0.1-8
#library(grid) # base
# Define a function to prepare silhouette images for plotting 
#prep.image <- function(string){
#	string <- readPNG(string)
#	string[,,2][which(string[,,2]==1)]<-NA
#	#string[,,2][which(string[,,2]==0)]<-1
#	string[,,2][which(string[,,2]<0.1)]<-1
#	string[,,2][which(string[,,2]!=1)]<-NA
#	string.raster <- as.raster(string[,,2],interpolate=F)
#	string.raster[string.raster=="#FFFFFF"] <- "#AAAAAA" # "#CCCCCC"
#	string.g <- rasterGrob(string.raster,interpolate=T)
#	return(string.g)
#}

# Only use the above function for adding silhouettes (not provided in public deposition)

#Balearica.g <- prep.image('Balearica_regulorum_silhouette.png')
#Upupa.g <- prep.image('Upupa_epops_silhouette.png')
#Topaza.g <- prep.image('Topaza_pyra_silhouette.png')
#Pelecanus.g <- prep.image('Pelecanus_occidentalis_silhouette.png')
#Pelagodroma.g <- prep.image('Pelagodroma_marina_silhouette.png')
#Nyctibius.g <- prep.image('Nyctibius_griseus_silhouette.png')
#Leucocarbo.g <- prep.image('Leucocarbo_atriceps_silhouette.png')
#Jynx.g <- prep.image('Jynx_torquilla_silhouette.png')
#Fregata.g <- prep.image('Fregata_silhouette.png')
#Columba.g <- prep.image('Columba_palumbus_silhouette.png')
# Load and prepare images

# Subclade ternary diagrams:

Cursorimorphae<-
ggplot( ) +
geom_segment(data=lines,aes(x=X1,xend=X3,y=X2,yend=X4),linewidth=fineline.width,col='grey',alpha=1)+
geom_path(data=outline,aes(x=x,y=y),linewidth=fineline.width,col='black')+
geom_segment(data=ticks,aes(x=X1,xend=X3,y=X2,yend=X4),linewidth=fineline.width,col='black',alpha=1)+
geom_path(data=dfline,aes(x=x,y=y),size=fineline.width*2,color='black',arrow = arrow(angle = 15, ends = "both", type = "closed",length=unit(0.2,'cm')))+
geom_path(data=dfline2,aes(x=x,y=y),size=fineline.width*2,color='black',arrow = arrow(angle = 15, ends = "both", type = "closed",length=unit(0.2,'cm')))+
#annotation_custom(grob=Balearica.g, xmin=-58, xmax=-33, ymin=32,ymax=45)+
geom_point(data=dfh[which(df$clade=='Cursorimorphae'),],aes(x=x,y=y),size=point.size/1.5,color='black')+
labs(fill = 'Development'    )+
geom_point(data=dfh[which(df$clade=='Cursorimorphae'),][rev(order(dfh[which(df$clade=='Cursorimorphae'),]$mode)),],aes(x=x,y=y, fill=factor(mode) ),size=fill.size/1.5,color='NA',shape=21)+
scale_discrete_manual(breaks=unique(df$mode)[c(8,4,7,3,6,2,5,1)], aesthetics='fill',values=c(rev(rainbow(100)[c(1,10,15,65,60,58,56,50)]))[c(5,1,4,6,3,7,2,8)] )+
geom_text(aes(x= -(32.5 + 35*cos(deg2rad(60))), y= 35*sin(deg2rad(60)) ) , label='PC1',angle=330, size=fontsize)+
geom_text(aes(x= -(32.5 + 46*cos(deg2rad(60))), y= 46*sin(deg2rad(60)) ) , label='PC2',angle=60,color='black', size=fontsize)+
scale_x_continuous(limits=range(c(ticks[,1],ticks[,3])),expand=c(0,0))+
theme(  panel.background = element_rect(fill='transparent'),
      plot.background = element_rect(fill='transparent', color=NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
	axis.title.x=element_blank(),
	axis.title.y=element_blank(),
	axis.text.x=element_blank(),
	axis.text.y=element_blank(),
	axis.ticks.y=element_blank(),
	axis.ticks.x=element_blank(),
	legend.position='none' )+
 coord_fixed()


Phaethoquornithes<-
ggplot( ) +
geom_segment(data=lines,aes(x=X1,xend=X3,y=X2,yend=X4),linewidth=fineline.width,col='grey',alpha=1)+
geom_path(data=outline,aes(x=x,y=y),linewidth=fineline.width,col='black')+
geom_segment(data=ticks,aes(x=X1,xend=X3,y=X2,yend=X4),linewidth=fineline.width,col='black',alpha=1)+
geom_path(data=dfline,aes(x=x,y=y),size=fineline.width*2,color='black',arrow = arrow(angle = 15, ends = "both", type = "closed",length=unit(0.2,'cm')))+
geom_path(data=dfline2,aes(x=x,y=y),size=fineline.width*2,color='black',arrow = arrow(angle = 15, ends = "both", type = "closed",length=unit(0.2,'cm')))+
#annotation_custom(grob=Pelagodroma.g, xmin=-54, xmax=-33, ymin=30,ymax=45)+
geom_point(data=dfh[which(df$clade=='Phaethoquornithes'),],aes(x=x,y=y),size=point.size/1.5,color='black')+
labs(fill = 'Development'    )+
geom_point(data=dfh[which(df$clade=='Phaethoquornithes'),][rev(order(dfh[which(df$clade=='Phaethoquornithes'),]$mode)),],aes(x=x,y=y, fill=factor(mode) ),size=fill.size/1.5,color='NA',shape=21)+
scale_discrete_manual(breaks=unique(df$mode)[c(8,4,7,3,6,2,5,1)], aesthetics='fill',values=c(rev(rainbow(100)[c(1,10,15,65,60,58,56,50)]))[c(5,1,4,6,3,7,2,8)] )+
geom_text(aes(x= -(32.5 + 35*cos(deg2rad(60))), y= 35*sin(deg2rad(60)) ) , label='PC1',angle=330, size=fontsize)+
geom_text(aes(x= -(32.5 + 46*cos(deg2rad(60))), y= 46*sin(deg2rad(60)) ) , label='PC2',angle=60,color='black', size=fontsize)+
scale_x_continuous(limits=range(c(ticks[,1],ticks[,3])),expand=c(0,0))+
theme(  panel.background = element_rect(fill='transparent'),
      plot.background = element_rect(fill='transparent', color=NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
	axis.title.x=element_blank(),
	axis.title.y=element_blank(),
	axis.text.x=element_blank(),
	axis.text.y=element_blank(),
	axis.ticks.y=element_blank(),
	axis.ticks.x=element_blank(),
	legend.position='none' )+
 coord_fixed()


Telluraves<-
ggplot( ) +
geom_segment(data=lines,aes(x=X1,xend=X3,y=X2,yend=X4),linewidth=fineline.width,col='grey',alpha=1)+
geom_path(data=outline,aes(x=x,y=y),linewidth=fineline.width,col='black')+
geom_segment(data=ticks,aes(x=X1,xend=X3,y=X2,yend=X4),linewidth=fineline.width,col='black',alpha=1)+
geom_path(data=dfline,aes(x=x,y=y),size=fineline.width*2,color='black',arrow = arrow(angle = 15, ends = "both", type = "closed",length=unit(0.2,'cm')))+
geom_path(data=dfline2,aes(x=x,y=y),size=fineline.width*2,color='black',arrow = arrow(angle = 15, ends = "both", type = "closed",length=unit(0.2,'cm')))+
#annotation_custom(grob=Upupa.g, xmin=-50, xmax=-30, ymin=26,ymax=45)+
geom_point(data=dfh[which(df$clade=='Telluraves'),],aes(x=x,y=y),size=point.size/1.5,color='black')+
labs(fill = 'Development'    )+
geom_point(data=dfh[which(df$clade=='Telluraves'),][rev(order(dfh[which(df$clade=='Telluraves'),]$mode)),],aes(x=x,y=y, fill=factor(mode) ),size=fill.size/1.5,color='NA',shape=21)+
scale_discrete_manual(breaks=unique(df$mode)[c(8,4,7,3,6,2,5,1)], aesthetics='fill',values=c(rev(rainbow(100)[c(1,10,15,65,60,58,56,50)]))[c(5,1,4,6,3,7,2,8)] )+
geom_text(aes(x= -(32.5 + 35*cos(deg2rad(60))), y= 35*sin(deg2rad(60)) ) , label='PC1',angle=330, size=fontsize)+
geom_text(aes(x= -(32.5 + 46*cos(deg2rad(60))), y= 46*sin(deg2rad(60)) ) , label='PC2',angle=60,color='black', size=fontsize)+
scale_x_continuous(limits=range(c(ticks[,1],ticks[,3])),expand=c(0,0))+
theme(  panel.background = element_rect(fill='transparent'),
      plot.background = element_rect(fill='transparent', color=NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
	axis.title.x=element_blank(),
	axis.title.y=element_blank(),
	axis.text.x=element_blank(),
	axis.text.y=element_blank(),
	axis.ticks.y=element_blank(),
	axis.ticks.x=element_blank(),
	legend.position='none' )+
 coord_fixed()


Strisores<-
ggplot( ) +
geom_segment(data=lines,aes(x=X1,xend=X3,y=X2,yend=X4),linewidth=fineline.width,col='grey',alpha=1)+
geom_path(data=outline,aes(x=x,y=y),linewidth=fineline.width,col='black')+
geom_segment(data=ticks,aes(x=X1,xend=X3,y=X2,yend=X4),linewidth=fineline.width,col='black',alpha=1)+
geom_path(data=dfline,aes(x=x,y=y),size=fineline.width*2,color='black',arrow = arrow(angle = 15, ends = "both", type = "closed",length=unit(0.2,'cm')))+
geom_path(data=dfline2,aes(x=x,y=y),size=fineline.width*2,color='black',arrow = arrow(angle = 15, ends = "both", type = "closed",length=unit(0.2,'cm')))+
#annotation_custom(grob=Topaza.g, xmin=-54, xmax=-32, ymin=26,ymax=46)+
geom_point(data=dfh[which(df$clade=='Strisores'),],aes(x=x,y=y),size=2*point.size/1.5,color='black')+
labs(fill = 'Development'    )+
geom_point(data=dfh[which(df$clade=='Strisores'),][rev(order(dfh[which(df$clade=='Strisores'),]$mode)),],aes(x=x,y=y, fill=factor(mode) ),size=2*fill.size/1.5,color='NA',shape=21)+
scale_discrete_manual(breaks=unique(df$mode)[c(8,4,7,3,6,2,5,1)], aesthetics='fill',values=c(rev(rainbow(100)[c(1,10,15,65,60,58,56,50)]))[c(5,1,4,6,3,7,2,8)] )+
geom_text(aes(x= -(32.5 + 35*cos(deg2rad(60))), y= 35*sin(deg2rad(60)) ) , label='PC1',angle=330, size=fontsize)+
geom_text(aes(x= -(32.5 + 46*cos(deg2rad(60))), y= 46*sin(deg2rad(60)) ) , label='PC2',angle=60,color='black', size=fontsize)+
scale_x_continuous(limits=range(c(ticks[,1],ticks[,3])),expand=c(0,0))+
theme(  panel.background = element_rect(fill='transparent'),
      plot.background = element_rect(fill='transparent', color=NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
	axis.title.x=element_blank(),
	axis.title.y=element_blank(),
	axis.text.x=element_blank(),
	axis.text.y=element_blank(),
	axis.ticks.y=element_blank(),
	axis.ticks.x=element_blank(),
	legend.position='none' )+
 coord_fixed()


library(ggpubr) # 0.6.0
left <- ggarrange(subplota,subplotb,ncol=1,nrow=2,align='hv',labels=c('(a)','(b)'), vjust=2, font.label = list(size = 11, color = "black", face = "italic"))
clades <- ggarrange(Cursorimorphae,Phaethoquornithes,Telluraves,Strisores,ncol=1,align='hv',labels=c('Cursorimorphae','Phaethoquornithes',
'Telluraves','Strisores'), vjust=1.2, hjust=c(-.65,-.5,-1,-1.25), font.label = list(size = 11, color = "black", face = "italic")) 
final <-ggarrange(left,clades,ncol=2,widths=c(2/3,1/3),labels=c('','(c)'), vjust=3, font.label = list(size = 11, color = "black", face = "italic") )

ggsave(filename='Figure_3_07_31_2025.pdf', final, width=18,height=14,unit='cm',device="pdf",family='Times') 
# Save the plot


# Script concludes 
