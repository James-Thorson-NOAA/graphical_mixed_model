
library(RTMB)
library(Matrix)
library(mvtnorm)
library(ape)
library(sem)
library(ggplot2)
library(gridExtra)

#root_dir = R'(C:\Users\James.Thorson\Desktop\Git\GGMM)'
root_dir = R'(C:\Users\James.Thorson\Desktop\Git\graphical_mixed_model)'
data_dir = file.path( root_dir, "data" )

# Extract single tree from VertTree, to save in github
if( FALSE ){
  local_dir = R'(C:\Users\James.Thorson\Desktop\Work files\AFSC\2025-04 -- Mammal tree)'
  tree = read.nexus( file.path( local_dir, "output.nex" ) )
  tree = tree$tree_1998
  write.tree( tree,
              file = file.path(data_dir,"VertTree_mammals.tre") )
}

#
Date = Sys.Date()
run_dir = file.path(root_dir, Date )
  dir.create(run_dir)

# Functions
source( file.path(root_dir, "trait_functions.R") )

# Load data
tree = read.tree( file.path( data_dir, "VertTree_mammals.tre" ) )
max_edge = max(tree$edge.length)
tree$edge.length = tree$edge.length / max_edge
s_root = Ntip(tree) + 1
n_nodes = Nnode(tree)
n_tips = Ntip(tree)

traits = read.table( file.path(data_dir,"PanTHERIA_1-0_WR05_Aug2008.txt"),
            row.names = NULL )

tree$tip.label = tolower( tree$tip.label )
traits$traits_binom = paste0( tolower(traits$MSW05_Species), "_", traits$MSW05_Binomial )

trait_to_tip = match( traits$traits_binom, tree$tip.label )
traits = traits[ which(!is.na(trait_to_tip)), ]

add_NA = function(vec) ifelse( as.numeric(vec) == -999, NA, as.numeric(vec) )
data = data.frame(
  ln_metabolism = log(add_NA(traits[,'X18.1_BasalMetRate_mLO2hr'])),
  ln_range = log(add_NA(traits[,'X22.1_HomeRange_km2'])),
  ln_size = log(add_NA(traits[,'X5.1_AdultBodyMass_g']))
)
rownames(data) = traits$traits_binom

#
sem = "
  ln_size -> ln_metabolism, b1
  ln_size -> ln_range, b2
"
ou_j = c(TRUE, TRUE, TRUE)

#
plm = phylolm::phylolm(
  ln_metabolism ~ ln_size,
  data = data,
  phy = tree,
  model = "OUrandomRoot"
)

# Table of data
has_data = ifelse( is.na(data), 0, 1 )
count = t(has_data) %*% has_data
write.csv( count, file = file.path(run_dir,"traits_data_available.csv") )


#############
# Fit in RTMB
#############

# Whether to drop root for ISAR
drop_bm_root = TRUE
which_drop = ((n_tips+n_nodes) * (seq_along(ou_j)-1) + (n_tips+1))[(ou_j==FALSE)]
assemble_version = 3

#
if( is.null(tree$node.label) & (n_nodes > 0) ){
  tree$node.label = paste0("node_",seq_len(n_nodes))
}

#
SEM_model = specifyModel( text=sem,
                          exog.variances=TRUE,
                          endog.variances=TRUE,
                          covs=colnames(data),
                          quiet=TRUE )
model = build_ram( model = SEM_model,
           vars = colnames(data) )

y_sj = as.matrix(data[match(c(tree$tip.label,tree$node.label),rownames(data)),,drop=FALSE])
parlist = list(
  y_sj = y_sj,
  beta_p = rep(1,max(model$parameter)),
  ln_theta = rep(0,ncol(data)),
  xbar = rep(0,ncol(data))
)
map = list(
  #ln_theta = factor( rep(1,ncol(data)) ),
  y_sj = ifelse( is.na(parlist$y_sj), seq_len(prod(dim(parlist$y_sj))), NA ),
  ln_theta = factor(ifelse(ou_j, seq_len(ncol(data)), NA))
)
if(isTRUE(drop_bm_root)) map$y_sj[which_drop] = NA
map$y_sj = factor(map$y_sj)

parlist$y_sj = ifelse( is.na(parlist$y_sj), 0, parlist$y_sj )

#method = "GMRF"
obj = MakeADFun( func = get_nll,
                  random = "y_sj",
                  map = map,
                  parameters = parlist,
                  silent = TRUE )
opt = nlminb( obj$par, obj$fn, obj$gr,
              control = list() )
rep = obj$report()
sdrep = sdreport(obj)

# Compare with phylosem ... logLik matches exactly when shared ln_theta
# map$ln_theta = factor(c(1,1,1))
psem = phylosem::phylosem(
  data = data,
  sem = sem,
  estimate_ou = all(ou_j),
  tree = tree
)


# Compare SDs
c( opt$par['beta_p'], psem$opt$par['beta_z'] )

#############
# Make figure
#############

theta = exp( obj$env$parList()$ln_theta ) / max_edge
ln_theta_se = as.list(sdrep, what = "Std. Error")$ln_theta
theta_lo = theta * exp(-1.96 * ln_theta_se)
theta_hi = theta * exp(1.96 * ln_theta_se)
max_depth = max(node.depth.edgelength( tree ))
xset = seq(0, max_depth * max_edge, length = 1000)
FUN = function(p) exp(-p*xset)
yset = sapply( theta, FUN = FUN )

library(viridisLite)
library(igraph)

DF = data.frame( from=model$first,
                 to=model$second,
                 label=obj$env$parList()$beta_p )
DF = subset( DF, DF$from != DF$to )
DF$label = round(DF$label, digits=3)
f = function(char_vec) gsub(x=char_vec, pattern="ln_", replace = "")
variables = f(colnames(data))
#variables = paste0( variables, "(t)" )
DF$from = f(DF$from)
DF$to = f(DF$to)


# Create and plotgraph
pg <- graph_from_data_frame( d = DF,
                             directed = TRUE,
                             vertices = data.frame(variables) )

png( file = file.path(run_dir,"traits.png"), width=4, height = 5, res = 200, units = "in")
  par( mfrow = c(2,1), mar=c(0,2.5,1,0.5), mgp = c(1.5,0.25,0), tck = -0.02, oma = c(3,0,0,0) )

  # Panel-1
  coords = layout_(pg, with_sugiyama())
  plot( pg, layout = coords, vertex.color = "white", vertex.frame.color = NA,
        vertex.label = variables )
  title( "Trait interactions" )

  # Panel-2
  matplot(
    x = xset,
    y = yset,
    type = "l",
    lwd = 2,
    lty = "solid",
    col = viridisLite::viridis(3),
    xlab = "",
    ylab = "O-U correlation",
    xaxs = "i",
    yaxs = "i"
  )
  title( "Correlation over time" )
  mtext( side=1, text = "Evolutionary time (million years)", line=1.5 )
  legend( "topright", bty = "n", fill = viridis(3),
          legend = variables )
  for( i in seq_len(ncol(data)) ){
    polygon( x = c(xset, rev(xset)),
             y = c(FUN(theta_lo[i]), rev(FUN(theta_hi[i]))),
             col = viridis(3, alpha = 0.2)[i],
             border = NA )
  }
dev.off()

################################
# Simulation study
################################

library(phylosem)
library(phylolm)

rep = obj$report()
ln_theta0_r = obj$env$parList()$ln_theta
n_reps = 200
p_set = c( 0.5, 0.7, 0.9 )

rmse_rpjz = array(NA, dim = c(n_reps,length(p_set),3,2), 
                      dimnames = list("rep" = seq_len(n_reps), "missing" = p_set, "var" = colnames(y_sj), "model" = c("GGMM","PSEM")) )
lntheta_rpjmz = array(NA, dim = c(n_reps,length(p_set),3,2,2), 
                      dimnames = list("rep" = seq_len(n_reps), "missing" = p_set, "var" = colnames(y_sj), "model" = c("GGMM","PLM"), "value" = c("est","se")) )
for( r_i in seq_len(n_reps) ){
for( p_i in seq_along(p_set) ){
  set.seed(r_i + p_i * n_reps)
  p_missing = p_set[p_i] # proportion missing
  ysim_k = rmvnorm_prec( rep$Qprime_kk, n = 1, mu = rep$muprime_k )
  #yobs_sj = ifelse( is.na(y_sj), NA, array(ysim_k,dim=dim(y_sj)) )
  
  #
  yobs_sj = array(ysim_k, dim=dim(y_sj), dimnames=dimnames(y_sj))
  for(z in 1:3){
    missing1 = which(!(rownames(y_sj) %in% rownames(data)))
    missing2 = sample( which(rownames(y_sj) %in% rownames(data)), replace = FALSE, size=nrow(data)*p_missing )
    yobs_sj[union(missing1,missing2),z] = NA
  }
    
  #method = "GMRF"
  parlist$y_sj = ifelse( is.na(yobs_sj), 0, yobs_sj )
  map$y_sj = factor(ifelse( is.na(yobs_sj), seq_len(prod(dim(parlist$y_sj))), NA ))
  obj_r = RTMB::MakeADFun( func = get_nll,
                    random = "y_sj",
                    map = map,
                    parameters = parlist,
                    silent = TRUE,
                    DLL = "RTMB" )
  opt_r = nlminb( obj_r$par, obj_r$fn, obj_r$gr,
                control = list() )
  hat1_r = obj_r$env$parList()$y_sj[ match(rownames(data),rownames(y_sj)), ]
  sdrep_r = sdreport(obj_r)
  lntheta_rpjmz[r_i,p_i,,'GGMM','est'] = as.list(sdrep_r, what = "Estimate")$ln_theta
  lntheta_rpjmz[r_i,p_i,,'GGMM','se'] = as.list(sdrep_r, what = "Std. Error")$ln_theta
  
  # Compare with phylosem ... logLik matches exactly when shared ln_theta
  data_r = yobs_sj[ match(rownames(data),rownames(y_sj)), ]
  psem_r = phylosem(
    data = data_r,
    sem = sem,
    estimate_ou = all(ou_j),
    tree = tree,
    control = phylosem_control( getsd = FALSE )
  )
  hat2_r = psem_r$parhat$x_vj[ match(rownames(data),rownames(y_sj)), ]

  # phylolm in order of colnames(data)
  plm1 = phylolm::phylolm(
    ln_metabolism ~ ln_size,
    data = na.omit(as.data.frame(data_r[,c("ln_metabolism","ln_size"),drop=FALSE])),
    phy = tree,
    model = "OUrandomRoot"
  )
  plm2 = phylolm::phylolm(
    ln_range ~ ln_size,
    data = na.omit(as.data.frame(data_r[,c("ln_range","ln_size"),drop=FALSE])),
    phy = tree,
    model = "OUrandomRoot"
  )
  plm3 = phylolm::phylolm(
    ln_size ~ 1,
    data = na.omit(as.data.frame(data_r[,"ln_size",drop=FALSE])),
    phy = tree,
    model = "OUrandomRoot"
  )
  lntheta_rpjmz[r_i,p_i,,'phylolm','est'] = log(c( plm1$optpar, plm2$optpar, plm3$optpar ))

  #
  true_r = array(ysim_k,dim=dim(y_sj))[ match(rownames(data),rownames(y_sj)), ]
  rmse_rpjz[r_i,p_i,,1] = apply( hat1_r - true_r, FUN = \(v)sqrt(mean(v^2)), MARGIN = 2 )
  rmse_rpjz[r_i,p_i,,2] = apply( hat2_r - true_r, FUN = \(v)sqrt(mean(v^2)), MARGIN = 2 )
}}

# 
sim_results = list(
  rmse_rpjz = rmse_rpjz, 
  lntheta_rpjmz = lntheta_rpjmz
)
saveRDS( sim_results, file = file.path(run_dir,"sim_results.RDS") )

#
table1 = apply(rmse_rpjz, MARGIN = 2:4, FUN = mean, na.rm = TRUE)
df1 = expand.grid( dimnames(table1) )
df1$RMSE = as.vector(table1)
g1 = ggplot(df1) + 
  geom_point( aes(y = RMSE, x = missing, col = model), position=position_dodge(0.1) ) + 
  facet_grid( vars(var) )

#
error_rpjm = lntheta_rpjmz[,,,,'est'] - (rep(1,n_reps) %o% rep(1,length(p_set)) %o% ln_theta0_r %o% c(1,1))
table2 = apply( lntheta_rpjmz[,,,,'est'] - (rep(1,n_reps) %o% rep(1,length(p_set)) %o% ln_theta0_r %o% c(1,1)), 
                 MARGIN = 2:4, 
                 FUN = \(v)sqrt(mean(v^2,na.rm = TRUE)) )
df2 = expand.grid( dimnames(table2) )
df2$RMSE = as.vector(table2)
g2 = ggplot(df2) + 
  geom_point( aes(y = RMSE, x = missing, col = model), position=position_dodge(0.1) ) + 
  facet_grid( vars(var) )

# Plot together
df3 = rbind( cbind(df1,"value" = "trait"), cbind(df2,"value" = "log-rate") )
g3 = ggplot(df3) + 
  geom_bar( aes(y = RMSE, x = missing, col = model, fill = model), stat = "identity", position = "dodge" ) +   # , position=position_dodge(0.1)
  facet_grid( rows = vars(var), cols = vars(value) )

#
z_rpj = error_rpjm[,,,'GGMM'] / lntheta_rpjmz[,,,'GGMM','se']
p_rpj = pnorm( z_rpj )
df4 = expand.grid( dimnames(p_rpj) )
df4$pvalue = as.vector(p_rpj)
df4$zvalue = as.vector(z_rpj)
g4 = ggplot(df4) + 
  geom_histogram( aes(x = pvalue, fill = missing), bins = 10 ) +  # , after_stat(count / sum(count))
  facet_grid( vars(var) )

g5 = grid.arrange( arrangeGrob(g3, g4, ncol=2) )
ggsave( g5, file = file.path(run_dir,"simulation.png"), width = 8, height = 5 )

# Plot together
g6 = ggplot(df3) + 
  geom_bar( aes(y = RMSE, x = missing, col = model, fill = model), stat = "identity", position = "dodge" ) +   # , position=position_dodge(0.1)
  facet_grid( cols = vars(var), rows = vars(value), scale = "free_y" )
ggsave( g6, file = file.path(run_dir,"errors.png"), width = 6, height = 3 )

#
df5 = expand.grid( dimnames(rmse_rpjz) )
  df5$Error = as.vector(rmse_rpjz)
df6 = expand.grid( dimnames(error_rpjm) )
  df6$Error = as.vector(error_rpjm)
df7 = rbind( cbind(df5,"value" = "Trait"), cbind(df6,"value" = "log(Rate)") )
g7 = ggplot(df7) + 
  geom_violin( aes(x=missing, y = Error, fill = model, col = model) ) + 
  facet_grid( rows = vars(value), cols = vars(var), scale = "free_y" )+ 
  theme_minimal()
ggsave( g7, file = file.path(run_dir,"errors_violins.png"), width = 6, height = 3 )

#
#zmax = dnorm(0) * 0.5 * 200
z = seq(-4, 4, by = 0.01 )
df7 = expand.grid( list(z = z, "missing" = factor(p_set), "var" = colnames(y_sj)) )
df7$height = dnorm(z) * 0.5 * 200
g8 = ggplot() + 
  geom_histogram( data = df4, aes(x = zvalue, fill = missing), binwidth = 0.5 ) +  # , y = after_stat(count / sum(count))
  facet_grid( rows = vars(var), cols = vars(missing) ) + 
  geom_line( data = df7, aes(x = z, y = height) ) +
  labs(y = "Count of simulation replicates")
ggsave( g8, file = file.path(run_dir,"zvalue_histograms.png"), width = 6, height = 4 )
