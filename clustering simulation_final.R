library(MASS)
library(abind)
library(GUniFrac)
library(MiRKAT)
library(MiSPU)
library(dirmult)
library(MCMCpack)
library(vegan)
library(fBasics)
library(PearsonDS)


set.seed(1240)
data("throat.otu.tab")
data("throat.tree")
data("throat.meta")


diri_result<-dirmult(throat.otu.tab)
p.est = diri_result$pi
names(p.est) <- names(diri_result$pi)
theta <- diri_result$theta
otu.ids <- throat.tree$tip.label
p.est <- p.est[otu.ids]

#####cluster based on tree distance, the cluster selected as 14 and 16#####
num_cluster<-20
tree.dist <- cophenetic(throat.tree)
obj <- pam(tree.dist, num_cluster)
clustering <- obj$clustering
sum(p.est[which(clustering == 14)])
gplus <- (1 - theta) / theta
g.est <- p.est * gplus

#####parameter settings#####
n_sim<-1000
n_sample<-100
depth<-1000
a<-0
b<-0

#####some functions######
scale2 = function(x)as.numeric(scale(x))

#####formal simulation#####
p_MiRKAT<-rep(NA,n_sim)
p_MiSPU<-rep(NA,n_sim)
p_adaptive<-rep(NA,n_sim)
p_anova_Du<-rep(NA,n_sim)
p_anova_DBC<-rep(NA,n_sim)
p_anova_Dw<-rep(NA,n_sim)
p_multiple_anova<-rep(NA,n_sim)

comm <- matrix(0, n_sample, length(p.est))
rownames(comm) <- 1:nrow(comm)
colnames(comm) <- names(g.est)
comm.p <- comm
nSeq <- rnbinom(n_sample , mu = depth, size = 25)

####preparing OTU#####
for (i in 1:n_sample ) {
  comm.p[i, ] <- rdirichlet(1, g.est)[1, ]
  comm[i, ] <- rmultinom(1, nSeq[i], prob=comm.p[i, ])[, 1]
}

for (j in 1:n_sim) {
  ####preparing response####
  otu.ids <- names(which(clustering == 14))
  OTU = comm[, otu.ids]
  X   = cbind(as.numeric(runif(n_sample)<0.5), rnorm(n_sample) + scale2(rowSums(OTU))*a)
  y <- scale(rowSums(OTU))[, 1] * b + rnorm(n_sample) + 0.5*(X[,1] + X[,2])
  
  ####calculate distance####
  unifracs = GUniFrac::GUniFrac(comm, throat.tree, alpha = c(0,0.5,1))$unifracs
  Dw  = unifracs[,,"d_1"]
  Du  = unifracs[,,"d_UW"]
  D.BC= as.matrix(vegdist(comm/nSeq, method="bray"))
  dis_MiRKAT<-list(Dw, Du, D.BC)
  Ks = lapply(dis_MiRKAT, FUN = function(d) D2K(d))
  p_MiRKAT[j]<-MiRKAT(y, X, Ks = Ks, out_type="C", method = "davies")$omnibus_p
  p_MiSPU[j]<-MiSPU(y, comm, throat.tree, cov = X, model ="gaussian",
        pow = c(2:8, Inf), n.perm = 5000)$aMiSPU$pvalue
  dis_adaptive<-abind(Dw,D.BC,along = 3)
  p_adaptive[j]<-adaptive_PREANOVA(dis_adaptive, y, confonding_stat = T, 
                                         confonding_var = X, r_vec = c(0.125,0.25,0.5,1,2))$final_p_value
  p_anova_Dw[j]<-tradition_preanova(n_rep=1000,d_mat=Dw,predictor = y, confonding_stat = T,
                                    confonding_var = X)
  p_anova_DBC[j]<-tradition_preanova(n_rep=1000,d_mat=D.BC,predictor = y, confonding_stat = T,
                                    confonding_var = X)
  p_multiple_anova[j]<-multiple_distance_preanova(n_rep=1000,d_mat=dis_adaptive, predictor = y, confonding_stat = T,
                                    confonding_var = X)
  
  print(j)
}


sum(as.numeric(p_MiRKAT<0.05),na.rm = T)
sum(as.numeric(p_MiSPU<0.05),na.rm = T)
sum(as.numeric(p_adaptive<0.05),na.rm = T)
sum(as.numeric(p_anova_Dw<0.05),na.rm = T)
sum(as.numeric(p_anova_DBC<0.05),na.rm = T)
sum(as.numeric(p_multiple_anova<0.05),na.rm = T)


sum(as.numeric(p_MiRKAT<0.01),na.rm = T)
sum(as.numeric(p_MiSPU<0.01),na.rm = T)
sum(as.numeric(p_adaptive<0.01),na.rm = T)
sum(as.numeric(p_anova_Dw<0.01),na.rm = T)
sum(as.numeric(p_anova_DBC<0.01),na.rm = T)
sum(as.numeric(p_multiple_anova<0.01),na.rm = T)
