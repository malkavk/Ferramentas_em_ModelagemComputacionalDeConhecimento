
  set.seed(12345)
  n <- 10^2                           # Número de observações
  
  u <- 1                                # Aceleração idealizada para o experimento
  
  dt <- 0.1                             # Menor intervalo tempo medido no simulador
  
  measnoise <- 10                       # Erro de medida da posição
  accelnoise <- 3                   # Erro da aceleração
  
  simPearsonTes<-function(ob1,ob2){
    sn <- length(ob1)
    sum1 <- sum(ob1)
    sum2 <- sum(ob2)
    sum1sq <- sum(ob1^2)
    sum2sq <- sum(ob2^2)
    psum <- sum(ob1*ob2)
    num <- psum-(sum1*sum2)/sn
    den <- sqrt( (sum1sq - sum1^2/sn) * (sum2sq - sum2^2/sn))
    res <- 0
    if(den != 0){
      r <- num/den
    }
    return <- r
  }
  
  simEuclidianTes<-function(ob1,ob2){
    return <- 1/sum((ob1-ob2)^2)
  }
  
  compare <- function(){
    
    posthatXpost_DistanciaEuclidiana <- simEuclidianTes(poshat,pos)
    postmeasXpost_DistanciaEuclidiana <- simEuclidianTes(posmeas,pos)
    postmeasXposthat_DistanciaEuclidiana <- simEuclidianTes(posmeas,poshat)
    posthatXposthat_DistanciaPearson <- simPearsonTes(poshat,pos)
    postmeasXpos_DistanciaPearson <- simPearsonTes(posmeas,pos)
    postmeasXposthat_DistanciaPearson <- simPearsonTes(posmeas,poshat)
    
  }
  
  plotGraphs<-function(){
    
    aa <- (poshat-posmeas)^2/posmeas
    xx <- sum(aa)
    t <- 1:n
    
    par(mfcol=c(3,1))
    
    plot(posmeas,pch=".",col="RED")
    points(pos,pch=".")
    points(poshat,pch=".",col="GREEN",type="l")
    
    
    plot(t,aa,pch=".",type="l",col="RED");
    points(t,-sqrt(Sk)^2,pch=".",col="BLUE")
    points(t,sqrt(Sk)^2,pch=".",col="GREEN")
    legend("bottomright",legend=c('Erros de Estimação(Inovação)','Limite Superior do Intervalo de Confiança','Limite Inferior do Intervalo de Confiança'),,fil=c("RED","BLUE","GREEN"))
    
    plot(t, Innt^2,pch=".",type="l")
    
    
    
  }
  
  a <- matrix(data=c(1, 0, dt, 1), ncol=2,nrow=2)
  b <- matrix(data=c(dt^2/2, dt), nrow=2, ncol=1)
  c <- matrix(data=c(1,1),nrow=1,ncol=2)
  x <- matrix(data=c(0,0),nrow=2,ncol=1)# Posição inicial
  xhat <- matrix(data=c(14,14),nrow=2,ncol=1)                             # Posição estimada inicial
  Sz <- measnoise^2                     # Variância do ruído da medida
  Sw <- accelnoise^2*matrix(data=c(dt^4, dt^3/2, dt^3/2, dt^2), 
                            nrow=2, ncol=2) # Estimativa
  P <- matrix(data=c(15,15,15,15),nrow=2,ncol=2)     # Matriz de covariância de estado (pode ser com qualquer coisa)
  Pr <- array(dim=c(n,2,2))
  pos <- array(dim=n)                    # Posição idealizada (durante todo o experimento)
  pos[1] <-0
  poshat <- array(dim=n)                  # Posição estimada (durante todo o experimento)
  poshat[1] <- 0 
  posmeas <- array(dim=n)                # Posição medida (durante todo o experimento)
  posmeas[1] <- 0
  vel <- array(dim=n)                    # Velocidade idealizada (durante todo o experimento)
  vel[1] <- 0
  velhat <- array(dim=n)                 # Velocidade estimada (durante todo o experimento)
  velhat[1] <- 0
  zz <- array(dim=n)                     # Posição + ruído (durante todo o experimento)
  zz[1] <- 0
  
  Sk <- array(dim=n)
  qk <- array(dim=n)
  Innt <- array(dim=n)
  
  for(t in 1:n){
    ProcessNoise <- accelnoise * matrix(data=c(dt^2/2*rnorm(n=1), dt*rnorm(n=1)),
                                        nrow=2, ncol=1) # Ruído
    x <- a%*%x + b*u + ProcessNoise   # Posição "correta" + ruído
    MeasNoise <- measnoise * rnorm(n=1,sd = accelnoise) # Ruído da medição
    y <- c%*%x + MeasNoise              # Medição do ruído
    xhat <- a%*%xhat + b%*%u            # Estimação da posição
    Inn <- y - c%*%xhat                 # Matriz erro (innovation)
    
    s <- c%*%P%*%t(c) + Sz
    
    K <- a%*%P%*%t(c)%*%solve(s)        # Ganho
    
    xhat <- xhat + K%*%Inn              # Correção da estimação
    Pr[t,,] <- P
    P <- a%*%P%*%t(a) - a%*%P%*%t(c)%*%solve(s)%*%c%*%P%*%t(a) + 
      Sw # Matriz de covariância (grau de incerteza das estimativas)
    
    Sk[t] <- c%*%P%*%t(c)+Sz
    
    Innt[t] <- Inn
    
    qk[t] <- Inn%*%solve(c%*%P%*%t(c) +Sz)%*%Inn
    
    pos[t] <- x[1]
    posmeas[t] <- y
    poshat[t] <- xhat[1]
    vel[t] <- x[2]
    velhat[t] <- xhat[2]
    
    z <- a%*%x + b%*%u # Posição sem ruído (idealizada)
    zz[t] <- z[1]
  }
  
  plotGraphs()
  
  
  poshatComp<-chisq.test(matrix(data=c(pos-min(posmeas),poshat-min(posmeas)),nrow=n,ncol=2))
  
  
  posmeasComp<-chisq.test(matrix(data=c(posmeas-min(posmeas),pos-min(posmeas)),nrow=n,ncol=2))
 dim(posmeas) 
 length(posmeas)