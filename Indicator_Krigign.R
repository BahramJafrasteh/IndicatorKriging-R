# Indicator kriging using gstat output

# Bahram Jafrasteh
# Ph.D. candidate
# Mining Engineering Department
# Isfahan University of Technology
# Isfahan, Iran.

# October 2017

IndicatorKriging <- function(DataFileName, Indexes_Inputs, 
                             CDFThs, DataFormat, OutFileName)
{
  
  # train data
  if ( isEqual(DataFormat, "Matlab") )
  {
    DF <- ReadMatlab(DataFileName, Indexes_Inputs)
  }
  else if (isEqual(DataFormat, "csv") )
  {
    DF <- ReadCSV(DataFileName, Indexes_Inputs)
  }
  else
  {
    stop("Unsupported Data Foramt. \n
         Data Format should be csv or Matlab. \n")
  }
  list2env( DF,env=environment())
  TrainD <- data.frame(Xtr,Ytr)
  # test data
  TestD <- data.frame(Xt,Yt)
  Dimension = dim(Xtr)[2]
  
  if (Dimension == 2)
  {
    coordinates(TrainD) <- ~X1+X2
    coordinates(TestD) <- ~X1+X2      
  }
  else if (Dimension == 3)
  {
    coordinates(TrainD) <- ~X1+X2+X3
    coordinates(TestD) <- ~X1+X2+X3
  }
  
  projection(TrainD) <- CRS("+init=epsg:3395")
  projection(TestD) <- CRS("+init=epsg:3395")
  
  # variogram
  CDFv = ecdf(TrainD$Ytr)
  Sequneces = seq(0, max(TrainD$Ytr), by = 0.01)
  CDFVals = CDFv(Sequneces)
  print(plot(CDFv))
  
  Len_Ths = length(CDFThs)
  
  IndicatorVals = matrix(list(), Len_Ths,1)
  
  ValThs=matrix(NA, nrow = Len_Ths-1 , ncol = 1)
  
  IndicatorVals$K =  matrix(NA,nrow = length(Yt), ncol = Len_Ths)
  
  IndicatorVals$Kv =  matrix(NA,nrow = length(Yt), ncol = Len_Ths)
  
  for (m in 1:Len_Ths){
    dis =  abs(CDFVals - CDFThs[m])
    ind = which.min(dis)
    ValThs [m] = Sequneces[ind]
    cat("The cCDFv is ", CDFThs[m], " and the corresponding value is ",
        ValThs[m], "\n")
    TrainD$indK <- (TrainD$Ytr <= ValThs[m])
    print(summary(TrainD$indK))
    
    
    KeepVarPars = "N"
    isSatisfied = "N"
    VarNug = 0.02
    VarRange = 200
    VarModelName = "Sph"
    VarSill = 0.2
    if (Dimension == 2)
    {
      # principal direction of anisotropy (clockwise from y) |
      A_p = 0
      # anisotropy ratio
      A_rat = 1
      Anisotropy = c(A_p, A_rat) 
    }
    else if (Dimension == 3)
    {
      # p: angle of the principal direction| 
      A_p = 0
      # q: Dip Angle of the principal Direction|
      A_q = 0
      # r: Third rotation angle|
      A_r = 0
      # s: major range/minor range|
      A_s = 1
      # t: major range/minor range|
      A_t = 1
      
      Anisotropy = c(A_p,A_q,A_r,A_s,A_t)
    }
    
    if (CDFThs[m] >= 0.9 ){LagDist = 10
    cutoffDist = 150;
    }
    else {LagDist = 20
    cutoffDist = 400}
    while ( isEqual (isSatisfied, "N") )
    {
      
      
      while ( isEqual (KeepVarPars, "N") )
      {
        if (Dimension == 2)
        {
          Variogm <- variogram(indK ~ 1, location = TrainD
                               , LagDistth = LagDist, cutoff = cutoffDist,
                               alpha = c(A_p, A_p+45, A_p+90, A_p+135) ,
                               beta = 0)
          
          
        }
        else if (Dimension == 3)
        {
        Variogm <- variogram(indK ~ 1, location = TrainD
                             , LagDistth = LagDist, cutoff = cutoffDist,
                             alpha = c(A_p, A_p+45, A_p+90, A_p+135) ,
                             beta = A_q
                             )
        }
        print(plot(Variogm, main =
                     paste("Variogram for threshold",ValThs[m], "CDF value", CDFThs[m] ),
                   ylab = 'gamma',xlab = 'lag distance', col ="black", lwd = 2) )
        ac = readChar(
          "Do you accept the plotted variogram[ Y:(Yes), N:(No) ]?",
          c("Y", "N"))
        if (isEqual( ac, "Y")) {
          KeepVarPars = "Y"}
        else{
          cat("The old parameters are: \n",
              "lag distance: ", LagDist, "\n The largerst distance: ", cutoffDist,  " \n")
          
          LagDist = readdouble(
            "Please enter new lag distance : ")
          cutoffDist = readdouble(
            "Please enter the largest distance for the plot : ")
        }
      }
      
      
      while (isEqual( isSatisfied, "N")){
        FittedVar <- fit.variogram(Variogm, vgm(VarSill, VarModelName, VarRange
                                                ,VarNug, ains = Anisotropy))
        if (Dimension == 2)
        {
          cat("Anisotropy parameters :",
              "\n \t Angle of the principal direction : ", A_p, 
              "\n \t Anisotropy ration :", A_rat, "\n")}
        else if (Dimension == 3)
        {
          cat("Anisotropy parameters :",
              "\n \t Angle of the principal direction : ", A_p, 
              "\n \t Dip Angle of the principal Direction :", A_q,
              "\n \t Third Rotation Angle : ", A_r,
              "\n \t Major range/Minor range : ", A_s,
              "\n \t Major range/Minor range : ", A_t, "\n")
        }
        
        print(plot(Variogm, model = FittedVar, main =
                     paste("Variogram model for threshold",ValThs[m], "CDF value", CDFThs[m] ),
                   ylab = 'gamma',xlab = 'lag distance', col ="black", lwd = 2))
        print(FittedVar)
        
        if (FittedVar$psill[2] == 0 ){
          isSatisfied = "N" 
        }
        else{ 
          #isSatisfied = 1}
          isSatisfied =  readChar("Are you satisfied by the provided model?[ Y:(Yes), N:(No) ]? ",
                                  c("Y", "N"))
        }
        if (isEqual( isSatisfied, "N")){
          if (Dimension == 2)
          {
            cat("The old parameters of variogram model are as follows \n",
                "Model Name: ", VarModelName, "\n Model Nugget: ", 
                VarNug,"\n Model Range: ", VarRange, " \n Model Sill: ",VarSill, " \n",
                "Anisotropy parameters :",
                "\n \t Angle of the principal direction : ", A_p, 
                "\n \t Anisotropy ration :", A_rat, "\n")  
          }
          else if (Dimension == 3){
            cat("The old parameters: \n",
                "Model Name: ", VarModelName, "\n Model Nugget: ", 
                VarNug,"\n Model Range: ", VarRange, " \n Model Sill: ",VarSill, " \n",
                "Anisotropy parameters :",
                "\n \t Angle of the principal direction : ", A_p, 
                "\n \t Dip Angle of the principal Direction :", A_q,
                "\n \t Third Rotation Angle : ", A_r,
                "\n \t Major range/Minor range : ", A_s,
                "\n \t Major range/Minor range : ", A_t, "\n")
          }
          cat("Previous Model Name: ", VarModelName, "\n")
          VarModelName = readChar("Please select a model name Sph:(Spherical), Exp:(Exponential) : \n",c("Sph", "Exp")
          )
          
          cat("Previous value of Nugget: ", VarNug, "\n")
          VarNug = readdouble("Please enter a new Nugget value : ")
          cat("Previous sill value : ", VarSill, "\n")
          VarSill = readdouble("Please enter a value for sill : ")
          cat("Previous range value : ", VarRange, "\n")
          VarRange = readdouble("Please enter a Range value : ")
          if (Dimension == 2)
          {
            cat("Previous Angle of the principal direction : ", A_p, "\n")
            A_p = readdouble("Angle of the principal direction : ")
            cat("Previous Anisotropy ratio : ", A_rat, "\n")
            A_rat = readdouble("Anisotropy ratio : ")
            Anisotropy = c(A_rng, A_rat)
          }
          else if (Dimension == 3)
          {
            cat("Previous Angle of the principal direction : ", A_p, "\n")
            A_p = readdouble("Angle of the principal direction : ")
            cat("Previous Dip Angle of the principal direction : ", A_q, "\n")
            A_q = readdouble("Dip Angle of the principal Direction : ")
            cat("Previous value of Third Rotation Angle : ", A_r, "\n")
            
            A_r = readdouble("Third Rotation Angle : ")
            cat("Major range/Minor range : ", A_s, "\n")
            A_s = readdouble("Major range/Minor range : ")
            cat("Major range/Minor range : ", A_t, "\n")
            A_t = readdouble("Major range/Minor range : ")
            Anisotropy = c(A_p, A_q, A_r, A_s, A_t)
          }
          
        }
        
        
        Sys.sleep(1)
        
      }
      if (FittedVar$range[2] > cutoffDist){
        isAllSet = 0; KeepVarPars = 0;
        isSatisfied = 0; 
        cutoffDist = cutoffDist/1.1; 
        LagDist = LagDist/1.1
      }
      else {
        KeepVarPars = "Y"
        isSatisfied = "Y"
      } 
      
    }
    
    # performing Ordinary Kriging of a threshold.
    TrainKrigTh <- krige(indK~1, TrainD,
                         TestD, model = FittedVar,
                         debug.level	= -1)
    # Recording estimated values with their variances
    IndicatorVals$K[,m] = TrainKrigTh$var1.pred
    IndicatorVals$Kv[,m] = TrainKrigTh$var1.var
    
  }
  # handling the order relation problem for IK.
  K = IndicatorVals$K
  for (ii in 1:length(Yt)){
    #ii = 1 
    for (jj in 1:Len_Ths){
      if (jj == 1) {
        K[ii,jj] = max(K[ii,jj],0)
      }
      else if (jj == Len_Ths){
        K[ii,jj] = max(K[ii,jj-1],K[ii,jj]);
        K[ii,jj] = min(K[ii,jj],1)  
      }
      
      else{
        K[ii,jj] = max(K[ii,jj-1],K[ii,jj])
      }
    }
    
    for (jj in Len_Ths:1){
      if (jj == 1) {
        K[ii,jj] = min(K[ii,jj],K[ii,jj+1])
        K[ii,jj] = min(K[ii,jj],0)
      }
      else if (jj == Len_Ths){
        K[ii,jj] = min(K[ii,jj],1)  
      }
      
      else{
        K[ii,jj] = min(K[ii,jj],K[ii,jj+1])
      }
    }
  }
  
  IndicatorVals$KT = cbind(matrix(0,nrow = length(Yt), ncol = 1),
                           K,matrix(1,nrow = length(Yt), ncol = 1))
  # Indicator Kriging E-type estimator.
  minx = min (TrainD$Ytr)
  maxx = max (TrainD$Ytr)
  ranges = cbind(minx,t(ValThs),maxx)
  # 
  Weights = matrix(1,nrow = Len_Ths+1, ncol = 1)
  # lower tail extrapolation power.
  Weights[1] = 1 #1.1
  # upper tail extrapolation power.
  Weights[Len_Ths+1] = 1 # 0.99
  
  IndicatorVals$KT.inv = matrix(NA,nrow = length(Yt), ncol = Len_Ths + 1) 
  # average of grade values.
  AveGrade =  matrix(NA,nrow = 1, ncol = Len_Ths + 1)
  
  for (l in 1:(Len_Ths+1)){
    
    Inds = which(Ytr <= ranges[l+1] & Ytr >= ranges[l] )
    
    AveGrade[l] = mean(Ytr[Inds])
    
    IndicatorVals$KT.inv[,l] = (IndicatorVals$KT[,l+1] - IndicatorVals$KT[,l])*
      (AveGrade[l])^Weights[l]
    
  }
  
  # E-type estimation
  IndicatorVals$IK = matrix(0,nrow = length(Yt), ncol = 1)
  
  
  for (j in 1:Len_Ths+1){
    IndicatorVals$IK = IndicatorVals$IK + IndicatorVals$KT.inv[,j]
    
  }
  # E-type conditional variance
  IndicatorVals$IKv = matrix(0,nrow = length(Yt), ncol = 1)
  
  for (kk in 2:Len_Ths + 2){
    temp = (AveGrade[kk-1] - IndicatorVals$IK)^2 * 
      (IndicatorVals$KT[,kk] - IndicatorVals$KT[,kk-1])
    IndicatorVals$IKv = IndicatorVals$IKv + temp
  }
  
  yHat_test = IndicatorVals$IK
  yHat_test_variance = IndicatorVals$IKv
  # calculate MSE of testing set
  ind_n = which(!is.na(yHat_test))
  MSE_test = mean( (yHat_test[ind_n]-Yt[ind_n])^2 )
  
  if ( !isEqual(OutFileName, "") )
  {
    cat ("Writing the Results to ", paste(OutFileName, ".csv", sep = "") )
    OutFileName = paste(OutFileName , ".csv", sep = "")
    write.csv(data.frame(yHat_test, yHat_test_variance), 
              file = OutFileName, append = FALSE, sep = ",",
              na = "NA", dec = ".", row.names = F,
              col.names = T)
  }
  
  
  
  list(yHat_test = yHat_test, yHat_test_variance = yHat_test_variance,
       MSE_test = MSE_test)
  
}

readdouble <- function(tex)
{ 
  n <- readline(prompt = paste(tex, "\t"))
  n <- as.double(n)
  if (is.na(n))
  {
    cat("Input ", n, " is not a double. Please type a double value.","\n")
    readdouble(n)
  }
  return(n)
}

readChar <- function(tex, states)
{ 
  
  n <- readline(prompt = paste(tex, "\t"))
  if (!is.character(n))
  {
    cat("Input ", n, " is not a character. Please type a character.","\n")
    readChar(n, states)
  }
  n <- as.character(n)
  for (i in 1:length(states))
  {
    if (n == states[i] || n == tolower( states[i] ) ||
        n == toupper( states[i] ) )
    {
      return(n);
    }
  }
  cat("Please provide a valid character among [", states,  "] \n")
  readChar(n, states)
}






# comparing two characters.
isEqual <- function(tex, stat)
{ 
  if (tex == stat || tolower(tex) == toupper(stat)
      || toupper(tex) == toupper(stat) ||
      tolower(tex) == tolower(stat) ||
      toupper(tex) == tolower(stat) )
  {
    return(T)
  }
  else
    return(F)
}
# Reading data from a matlab file
ReadMatlab <- function(InputDataFileName, Indexes)
{
  # load the dataset1
  DataFile <- readMat.default(InputDataFileName)
  
  Att = attributes(DataFile)
  
  Att_ind = which(Att$names == "Xtr")
  di_tr=sapply(DataFile[Att_ind ],dim)
  Xtr <- matrix(unlist(DataFile[Att_ind]), 
                ncol = di_tr[2] , nrow = di_tr[1])
  
  Att_ind = which(Att$names == "Ytr")
  di=sapply(DataFile[Att_ind],dim)
  Ytr <- matrix(unlist(DataFile[Att_ind]), ncol = di[2], nrow = di[1])
  
  Att_ind = which(Att$names == "Xt")
  di_t=sapply(DataFile[Att_ind],dim)
  Xt <- matrix(unlist(DataFile[Att_ind]), ncol = di_t[2], nrow = di_t[1])
  
  Att_ind = which(Att$names == "Yt")
  di=sapply(DataFile[Att_ind],dim)
  Yt <- matrix(unlist(DataFile[Att_ind]), ncol = di[2], nrow = di[1])
  
  # Define inputs
  if ( length(Indexes) == 2 )
  {
    Xtr <- cbind(Xtr[,Indexes[1]],Xtr[,Indexes[2]])
    Xt <- cbind(Xt[,Indexes[1]],Xt[,Indexes[2]])
  }
  else if ( length(Indexes) == 3 )
  {
    Xtr <- cbind(Xtr[,Indexes[1]],Xtr[,Indexes[2]],Xtr[,Indexes[3]])
    Xt <- cbind(Xt[,Indexes[1]],Xt[,Indexes[2]],Xt[,Indexes[3]])
  }
  else
  {
    stop("Number of input variables should be 2 or 3.\n")
  }
  
  list(Xtr= Xtr, Ytr = Ytr, Xt= Xt, Yt = Yt)
}




ReadCSV <- function(InputTrainDataFileName, InputTestDataFileName)
{
  # load the dataset1
  TrainDataFile <- read.csv(InputTrainDataFileName,
                            sep = ",", header = F)
  
  TestDataFile <- read.csv(InputTestDataFileName,
                           sep = ",", header = F)
  # Define inputs
  if ( length(TrainDataFile) == 3 && length(TestDataFile) == 3)
  {
    Xtr <- cbind(TrainDataFile$V1, TrainDataFile$V2)
    Ytr <- TrainDataFile$V3
    Xt <- cbind(TestDataFile$V1, TestDataFile$V2)
    Yt <- TestDataFile$V3
  }
  else if ( length(TrainDataFile) == 4 && length(TestDataFile) == 4)
  {
    Xtr <- cbind(TrainDataFile$V1, TrainDataFile$V2, TrainDataFile$V3)
    Ytr <- TrainDataFile$V4
    Xt <- cbind(TestDataFile$V1, TestDataFile$V2, TestDataFile$V3)
    Yt <- TestDataFile$V4
  }
  else
  {
    stop("Dimension mismatch. \n")
  }
  
  list(Xtr= Xtr, Ytr = Ytr, Xt= Xt, Yt = Yt)
}