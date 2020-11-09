rm(list=ls())
library('doMC')
registerDoMC(16)
comb <- function(x, ...) {
  lapply(seq_along(x),
         function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
}

#parameters
mean_1<-210
mean_2<-210
sd_1<-23
sd_2<-36

c_1o<-5
c_2o<-7
c_1u<-6
c_2u<-6

unit_1a<-4
unit_1b<-7
unit_1c<-8
unit_2a<-6
unit_2b<-5
unit_2c<-8

time_a<-2200
time_b<-2500
time_c<-3500

rang<-seq(0,5,by=0.25)



#sp
sp.iter<-500000

sp.func<-function(par){
  x_1<-par[1]
  x_2<-par[2]
  cost_1<-sum(c_1o*(x_1-demand_1_sp)*((x_1-demand_1_sp)>0)+c_1u*(demand_1_sp-x_1)*((x_1-demand_1_sp)<0))/sp.iter
  cost_2<-sum(c_2o*(x_2-demand_2_sp)*((x_2-demand_2_sp)>0)+c_2u*(demand_2_sp-x_2)*((x_2-demand_2_sp)<0))/sp.iter
  return(cost_1+cost_2)
}

ui<-cbind(c(-unit_1a,-unit_1b,-unit_1c),c(-unit_2a,-unit_2b,-unit_2c))
ci<-c(-time_a,-time_b,-time_c)


##mean
sp.aver<-20

sp_record_1<-matrix(nrow = 21,ncol = sp.aver)
for(j in 1:sp.aver){
  sol_sp_mean_1<-c()

  for (i in rang){
  demand_1_sp<-rnorm(sp.iter,mean_1*(1+0.01*i),sd_1)
  demand_2_sp<-rnorm(sp.iter,mean_2,sd_2)
  opti<-constrOptim(c(0,0), sp.func, grad=NULL, ui=ui, ci=ci)
  sol_sp_mean_1<-c(sol_sp_mean_1,opti$value)
  }
  sp_record_1[,j]<-sol_sp_mean_1
}
sp_mean_1<-rowMeans(sp_record_1, na.rm = FALSE, dims = 1)



sp_record_2<-matrix(nrow = 21,ncol = sp.aver)
for(j in 1:sp.aver){
  sol_sp_mean_2<-c()
  
  for (i in rang){
    demand_1_sp<-rnorm(sp.iter,mean_1,sd_1)
    demand_2_sp<-rnorm(sp.iter,mean_2*(1+0.01*i),sd_2)
    opti<-constrOptim(c(0,0), sp.func, grad=NULL, ui=ui, ci=ci)
    sol_sp_mean_2<-c(sol_sp_mean_2,opti$value)
  }
  sp_record_2[,j]<-sol_sp_mean_2
}
sp_mean_2<-rowMeans(sp_record_2, na.rm = FALSE, dims = 1)



##var

sp_record_3<-matrix(nrow = 21,ncol = sp.aver)
for(j in 1:sp.aver){
  sol_sp_var_1<-c()
  
  for (i in rang){
    demand_1_sp<-rnorm(sp.iter,mean_1,sd_1*(1-0.01*i))
    demand_2_sp<-rnorm(sp.iter,mean_2,sd_2)
    opti<-constrOptim(c(0,0), sp.func, grad=NULL, ui=ui, ci=ci)
    sol_sp_var_1<-c(sol_sp_var_1,opti$value)
  }
  sp_record_3[,j]<-sol_sp_var_1
}
sp_var_1<-rowMeans(sp_record_3, na.rm = FALSE, dims = 1)



sp_record_4<-matrix(nrow = 21,ncol = sp.aver)
for(j in 1:sp.aver){
  sol_sp_var_2<-c()
  
  for (i in rang){
    demand_1_sp<-rnorm(sp.iter,mean_1,sd_1)
    demand_2_sp<-rnorm(sp.iter,mean_2,sd_2*(1-0.01*i))
    opti<-constrOptim(c(0,0), sp.func, grad=NULL, ui=ui, ci=ci)
    sol_sp_var_2<-c(sol_sp_var_2,opti$value)
  }
  sp_record_4[,j]<-sol_sp_var_2
}
sp_var_2<-rowMeans(sp_record_4, na.rm = FALSE, dims = 1)



#Simulation
sim.time<-10000


#sim_500
iter.500<-500

iter.500.func<-function(par){
  x_1<-par[1]
  x_2<-par[2]
  cost_1<-sum(c_1o*(x_1-demand_1_500)*((x_1-demand_1_500)>0)+c_1u*(demand_1_500-x_1)*((x_1-demand_1_500)<0))/iter.500
  cost_2<-sum(c_2o*(x_2-demand_2_500)*((x_2-demand_2_500)>0)+c_2u*(demand_2_500-x_2)*((x_2-demand_2_500)<0))/iter.500
  return(cost_1+cost_2)
}

##mean

sim.500.mean.1<-foreach (j=1:sim.time) %:% foreach (i=rang) %do%{
  demand_1_500<-rnorm(iter.500,mean_1*(1+0.01*i),sd_1)
  demand_2_500<-rnorm(iter.500,mean_2,sd_2)
  opti<-constrOptim(c(0,0), iter.500.func, grad=NULL, ui=ui, ci=ci)
  record_lp<-opti$value
}

sim.500.mean.2<-foreach (j=1:sim.time) %:% foreach (i=rang) %do%{
  demand_1_500<-rnorm(iter.500,mean_1,sd_1)
  demand_2_500<-rnorm(iter.500,mean_2*(1+0.01*i),sd_2)
  opti<-constrOptim(c(0,0), iter.500.func, grad=NULL, ui=ui, ci=ci)
  record_lp<-opti$value
}

##var
sim.500.var.1<-foreach (j=1:sim.time) %:% foreach (i=rang) %do%{
  demand_1_500<-rnorm(iter.500,mean_1,sd_1*(1-0.01*i))
  demand_2_500<-rnorm(iter.500,mean_2,sd_2)
  opti<-constrOptim(c(0,0), iter.500.func, grad=NULL, ui=ui, ci=ci)
  record_lp<-opti$value
}

sim.500.var.2<-foreach (j=1:sim.time) %:% foreach (i=rang) %do%{
  demand_1_500<-rnorm(iter.500,mean_1,sd_1)
  demand_2_500<-rnorm(iter.500,mean_2,sd_2*(1-0.01*i))
  opti<-constrOptim(c(0,0), iter.500.func, grad=NULL, ui=ui, ci=ci)
  record_lp<-opti$value
}

#sim_1000
iter.1000<-1000

iter.1000.func<-function(par){
  x_1<-par[1]
  x_2<-par[2]
  cost_1<-sum(c_1o*(x_1-demand_1_1000)*((x_1-demand_1_1000)>0)+c_1u*(demand_1_1000-x_1)*((x_1-demand_1_1000)<0))/iter.1000
  cost_2<-sum(c_2o*(x_2-demand_2_1000)*((x_2-demand_2_1000)>0)+c_2u*(demand_2_1000-x_2)*((x_2-demand_2_1000)<0))/iter.1000
  return(cost_1+cost_2)
}

##mean

sim.1000.mean.1<-foreach (j=1:sim.time) %:% foreach (i=rang) %do%{
  demand_1_1000<-rnorm(iter.1000,mean_1*(1+0.01*i),sd_1)
  demand_2_1000<-rnorm(iter.1000,mean_2,sd_2)
  opti<-constrOptim(c(0,0), iter.1000.func, grad=NULL, ui=ui, ci=ci)
  record_lp<-opti$value
}

sim.1000.mean.2<-foreach (j=1:sim.time) %:% foreach (i=rang) %do%{
  demand_1_1000<-rnorm(iter.1000,mean_1,sd_1)
  demand_2_1000<-rnorm(iter.1000,mean_2*(1+0.01*i),sd_2)
  opti<-constrOptim(c(0,0), iter.1000.func, grad=NULL, ui=ui, ci=ci)
  record_lp<-opti$value
}

##var
sim.1000.var.1<-foreach (j=1:sim.time) %:% foreach (i=rang) %do%{
  demand_1_1000<-rnorm(iter.1000,mean_1,sd_1*(1-0.01*i))
  demand_2_1000<-rnorm(iter.1000,mean_2,sd_2)
  opti<-constrOptim(c(0,0), iter.1000.func, grad=NULL, ui=ui, ci=ci)
  record_lp<-opti$value
}

sim.1000.var.2<-foreach (j=1:sim.time) %:% foreach (i=rang) %do%{
  demand_1_1000<-rnorm(iter.1000,mean_1,sd_1)
  demand_2_1000<-rnorm(iter.1000,mean_2,sd_2*(1-0.01*i))
  opti<-constrOptim(c(0,0), iter.1000.func, grad=NULL, ui=ui, ci=ci)
  record_lp<-opti$value
}



#sim_2000
iter.2000<-2000

iter.2000.func<-function(par){
  x_1<-par[1]
  x_2<-par[2]
  cost_1<-sum(c_1o*(x_1-demand_1_2000)*((x_1-demand_1_2000)>0)+c_1u*(demand_1_2000-x_1)*((x_1-demand_1_2000)<0))/iter.2000
  cost_2<-sum(c_2o*(x_2-demand_2_2000)*((x_2-demand_2_2000)>0)+c_2u*(demand_2_2000-x_2)*((x_2-demand_2_2000)<0))/iter.2000
  return(cost_1+cost_2)
}

##mean

sim.2000.mean.1<-foreach (j=1:sim.time) %:% foreach (i=rang) %do%{
  demand_1_2000<-rnorm(iter.2000,mean_1*(1+0.01*i),sd_1)
  demand_2_2000<-rnorm(iter.2000,mean_2,sd_2)
  opti<-constrOptim(c(0,0), iter.2000.func, grad=NULL, ui=ui, ci=ci)
  record_lp<-opti$value
}

sim.2000.mean.2<-foreach (j=1:sim.time) %:% foreach (i=rang) %do%{
  demand_1_2000<-rnorm(iter.2000,mean_1,sd_1)
  demand_2_2000<-rnorm(iter.2000,mean_2*(1+0.01*i),sd_2)
  opti<-constrOptim(c(0,0), iter.2000.func, grad=NULL, ui=ui, ci=ci)
  record_lp<-opti$value
}

##var
sim.2000.var.1<-foreach (j=1:sim.time) %:% foreach (i=rang) %do%{
  demand_1_2000<-rnorm(iter.2000,mean_1,sd_1*(1-0.01*i))
  demand_2_2000<-rnorm(iter.2000,mean_2,sd_2)
  opti<-constrOptim(c(0,0), iter.2000.func, grad=NULL, ui=ui, ci=ci)
  record_lp<-opti$value
}

sim.2000.var.2<-foreach (j=1:sim.time) %:% foreach (i=rang) %do%{
  demand_1_2000<-rnorm(iter.2000,mean_1,sd_1)
  demand_2_2000<-rnorm(iter.2000,mean_2,sd_2*(1-0.01*i))
  opti<-constrOptim(c(0,0), iter.2000.func, grad=NULL, ui=ui, ci=ci)
  record_lp<-opti$value
}

#SA

#sa_500
iter.500<-500

objective.in<-c(rep(c_1o,iter.500),rep(c_2o,iter.500),rep(c_1u,iter.500),rep(c_2u,iter.500),0,0)/iter.500
con_t1<-c(rep(0,4*iter.500),unit_1a,unit_2a)
con_t2<-c(rep(0,4*iter.500),unit_1b,unit_2b)
con_t3<-c(rep(0,4*iter.500),unit_1c,unit_2c)

mac<-matrix(0L,nrow=4*iter.500,ncol=4*iter.500)
diag(mac)<-1
mac<-cbind(mac,matrix(c(rep(-1,iter.500),rep(0,iter.500),rep(1,iter.500),rep(0,iter.500),
                        rep(0,iter.500),rep(-1,iter.500),rep(0,iter.500),rep(1,iter.500)),ncol=2))

const.mat<-rbind(con_t1,con_t2,con_t3,mac)
const.dir<-c(rep("<=",3),rep(">=",4*iter.500))

sa.500<-foreach (j=1:sim.time,.packages='lpSolve') %dopar%{
  demand_1_500<-rnorm(iter.500,mean_1,sd_1)
  demand_2_500<-rnorm(iter.500,mean_2,sd_2)
  const.rhs<-c(time_a,time_b,time_c,-demand_1_500,-demand_2_500,demand_1_500,demand_2_500)
  opti<-lp(direction="min", objective.in, const.mat,const.dir, const.rhs,compute.sens=TRUE)
  
  ##mean
  pie.mean.1<-sum(-opti$duals[4:(3+iter.500)],opti$duals[(4+2*iter.500):(3+3*iter.500)])
  sol_500_mean_1_SA<-pie.mean.1*mean_1*0.01*0.25
  
  pie.mean.2<-sum(-opti$duals[(4+iter.500):(3+2*iter.500)],opti$duals[(4+3*iter.500):(3+4*iter.500)])
  sol_500_mean_2_SA<-pie.mean.2*mean_2*0.01*0.25
  
  ##var
  pie.var.1<-c(-opti$duals[4:(3+iter.500)],opti$duals[(4+2*iter.500):(3+3*iter.500)])
  sol_500_var_1_SA<-sum(pie.var.1*(mean_1-demand_1_500))*0.01*0.25
  
  pie.var.2<-c(-opti$duals[(4+iter.500):(3+2*iter.500)],opti$duals[(4+3*iter.500):(3+4*iter.500)])
  sol_500_var_2_SA<-sum(pie.var.2*(mean_2-demand_2_500))*0.01*0.25
  
  ##record
  re<-c(sol_500_mean_1_SA,sol_500_mean_2_SA,sol_500_var_1_SA,sol_500_var_2_SA)
}

#sa_1000
iter.1000<-1000

objective.in<-c(rep(c_1o,iter.1000),rep(c_2o,iter.1000),rep(c_1u,iter.1000),rep(c_2u,iter.1000),0,0)/iter.1000
con_t1<-c(rep(0,4*iter.1000),unit_1a,unit_2a)
con_t2<-c(rep(0,4*iter.1000),unit_1b,unit_2b)
con_t3<-c(rep(0,4*iter.1000),unit_1c,unit_2c)

mac<-matrix(0L,nrow=4*iter.1000,ncol=4*iter.1000)
diag(mac)<-1
mac<-cbind(mac,matrix(c(rep(-1,iter.1000),rep(0,iter.1000),rep(1,iter.1000),rep(0,iter.1000),
                        rep(0,iter.1000),rep(-1,iter.1000),rep(0,iter.1000),rep(1,iter.1000)),ncol=2))

const.mat<-rbind(con_t1,con_t2,con_t3,mac)
const.dir<-c(rep("<=",3),rep(">=",4*iter.1000))

sa.1000<-foreach (j=1:sim.time,.packages='lpSolve') %dopar%{
  demand_1_1000<-rnorm(iter.1000,mean_1,sd_1)
  demand_2_1000<-rnorm(iter.1000,mean_2,sd_2)
  const.rhs<-c(time_a,time_b,time_c,-demand_1_1000,-demand_2_1000,demand_1_1000,demand_2_1000)
  opti<-lp(direction="min", objective.in, const.mat,const.dir, const.rhs,compute.sens=TRUE)
  
  ##mean
  pie.mean.1<-sum(-opti$duals[4:(3+iter.1000)],opti$duals[(4+2*iter.1000):(3+3*iter.1000)])
  sol_1000_mean_1_SA<-pie.mean.1*mean_1*0.01*0.25
  
  pie.mean.2<-sum(-opti$duals[(4+iter.1000):(3+2*iter.1000)],opti$duals[(4+3*iter.1000):(3+4*iter.1000)])
  sol_1000_mean_2_SA<-pie.mean.2*mean_2*0.01*0.25
  
  ##var
  pie.var.1<-c(-opti$duals[4:(3+iter.1000)],opti$duals[(4+2*iter.1000):(3+3*iter.1000)])
  sol_1000_var_1_SA<-sum(pie.var.1*(mean_1-demand_1_1000))*0.01*0.25
  
  pie.var.2<-c(-opti$duals[(4+iter.1000):(3+2*iter.1000)],opti$duals[(4+3*iter.1000):(3+4*iter.1000)])
  sol_1000_var_2_SA<-sum(pie.var.2*(mean_2-demand_2_1000))*0.01*0.25
  
  ##record
  re<-c(sol_1000_mean_1_SA,sol_1000_mean_2_SA,sol_1000_var_1_SA,sol_1000_var_2_SA)
}



#sa_2000
iter.2000<-2000

objective.in<-c(rep(c_1o,iter.2000),rep(c_2o,iter.2000),rep(c_1u,iter.2000),rep(c_2u,iter.2000),0,0)/iter.2000
con_t1<-c(rep(0,4*iter.2000),unit_1a,unit_2a)
con_t2<-c(rep(0,4*iter.2000),unit_1b,unit_2b)
con_t3<-c(rep(0,4*iter.2000),unit_1c,unit_2c)

mac<-matrix(0L,nrow=4*iter.2000,ncol=4*iter.2000)
diag(mac)<-1
mac<-cbind(mac,matrix(c(rep(-1,iter.2000),rep(0,iter.2000),rep(1,iter.2000),rep(0,iter.2000),
                        rep(0,iter.2000),rep(-1,iter.2000),rep(0,iter.2000),rep(1,iter.2000)),ncol=2))

const.mat<-rbind(con_t1,con_t2,con_t3,mac)
const.dir<-c(rep("<=",3),rep(">=",4*iter.2000))

sa.2000<-foreach (j=1:sim.time,.packages='lpSolve') %dopar%{
  demand_1_2000<-rnorm(iter.2000,mean_1,sd_1)
  demand_2_2000<-rnorm(iter.2000,mean_2,sd_2)
  const.rhs<-c(time_a,time_b,time_c,-demand_1_2000,-demand_2_2000,demand_1_2000,demand_2_2000)
  opti<-lp(direction="min", objective.in, const.mat,const.dir, const.rhs,compute.sens=TRUE)
  
  ##mean
  pie.mean.1<-sum(-opti$duals[4:(3+iter.2000)],opti$duals[(4+2*iter.2000):(3+3*iter.2000)])
  sol_2000_mean_1_SA<-pie.mean.1*mean_1*0.01*0.25
  
  pie.mean.2<-sum(-opti$duals[(4+iter.2000):(3+2*iter.2000)],opti$duals[(4+3*iter.2000):(3+4*iter.2000)])
  sol_2000_mean_2_SA<-pie.mean.2*mean_2*0.01*0.25
  
  ##var
  pie.var.1<-c(-opti$duals[4:(3+iter.2000)],opti$duals[(4+2*iter.2000):(3+3*iter.2000)])
  sol_2000_var_1_SA<-sum(pie.var.1*(mean_1-demand_1_2000))*0.01*0.25
  
  pie.var.2<-c(-opti$duals[(4+iter.2000):(3+2*iter.2000)],opti$duals[(4+3*iter.2000):(3+4*iter.2000)])
  sol_2000_var_2_SA<-sum(pie.var.2*(mean_2-demand_2_2000))*0.01*0.25
  
  ##record
  re<-c(sol_2000_mean_1_SA,sol_2000_mean_2_SA,sol_2000_var_1_SA,sol_2000_var_2_SA)
}



save(sp_mean_1,sp_mean_2,sp_var_1,sp_var_2
     ,sim.500.mean.1,sim.500.mean.2,sim.500.var.1,sim.500.var.2
     ,sim.1000.mean.1,sim.1000.mean.2,sim.1000.var.1,sim.1000.var.2
     ,sim.2000.mean.1,sim.2000.mean.2,sim.2000.var.1,sim.2000.var.2
     ,sa.500,sa.1000,sa.2000,file = 'simCon.Rdata')
