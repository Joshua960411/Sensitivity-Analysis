rm(list = ls())
load('simCon.Rdata')
sim.time<-10000
rang<-seq(0,5,by=0.25)

get.diffmac<-function(Listoflist){
  diffmac<-matrix(nrow = sim.time,ncol = length(rang)-1)
  for (i in 1:sim.time){
    for(j in 1:length(rang)-1){
      diffmac[i,j]<-unlist(Listoflist[[i]][j+1])-unlist(Listoflist[[i]][j])
    }
  }
  return(diffmac)
}

get.meanarray<-function(mac){
  meanarray<-colMeans(mac, na.rm = FALSE, dims = 1)
  return(meanarray)
}

get.direction<-function(mac){
  dirmac<-mac>=0
  return(dirmac)
}


diff_mean_1<-sp_mean_1[2:21]-sp_mean_1[1:20]
diff_mean_2<-sp_mean_2[2:21]-sp_mean_2[1:20]
diff_var_1<-sp_var_1[2:21]-sp_var_1[1:20]
diff_var_2<-sp_var_2[2:21]-sp_var_2[1:20]

#percentage
per.re1000.mean.1<-get.meanarray(get.diffmac(sim.1000.mean.1))/diff_mean_1
per.re1000.mean.2<-get.meanarray(get.diffmac(sim.1000.mean.2))/diff_mean_2
per.re1000.var.1<-get.meanarray(get.diffmac(sim.1000.var.1))/diff_var_1
per.re1000.var.2<-get.meanarray(get.diffmac(sim.1000.var.2))/diff_var_2



