#########################################################################################
#Data:2019/08/12
#Author；Wang
##references: I use the plotTAD function from https://github.com/menghaowei/ngstools/blob/master/README.md author's E-mail: meng_howard@126.com
###this function is to plot hic plot for inter or intra of chromosomes, and the data is h5 from NCBI
##h5 format from NCBI:
#group            name       otype  dclass     dim
#0     / balance_factors H5I_DATASET   FLOAT    6207
#1     /   bin_positions H5I_DATASET INTEGER  x 6207
#2     /   chr_bin_range H5I_DATASET INTEGER    x 25
#3     /            chrs H5I_DATASET  STRING      25
#4     /    interactions H5I_DATASET   FLOAT  x 6207

#to see the particular information of h5 
see_hic_h5<-function(file){
  if(require(rhdf5)){
    print("the summary of h5 is")
    print(paste("----------------------------------"))
    c<-h5ls(file)
    print(c)
    #interactions<-h5read(file,"interactions")
    #deal with the NaN is in plotTAD
    #read the bin positions and get the resolution
    bins<-h5read(file,"bin_positions")
    resolution<-bins[2,2]-bins[2,1]
    print("resolutions:")
    print(paste("this hic data's resolution is",resolution/1000,"kb"))
    print(paste("----------------------------------"))
    #read the chrs
    chrs<-h5read(file,"chrs")
    #read the  chr_bins_range<-h5read
    chr_bins_range<-h5read(file,"chr_bin_range")
    colnames(chr_bins_range)<-chrs
    #print(paste("----------------------------------"))
    print("the chr bin range and chr range:")
    print("the chr bin range is ")
    print(chr_bins_range)
    chr_bins_range<-as.numeric( chr_bins_range)
    print("the chr range is ")
    print(chr_bins_range*resolution)
    print(paste("----------------------------------"))
  }
  else{
    print("please install rhdf5")
  }
}

#plot intra-chromosome interactions
ploth5_intra_hic<-function(file,chr,from,to,col.min = "red",col.max = "red",col.boundary = 0,ylim=NULL,minBound=0,maxBound=0.95,mat.upper=T,mat.part=T,title=TRUE){
  
  
  #file :h5 file from NCBI format is above
  #chr plot which chromosome
  #from to position such as chr1,1000000,1500000
  if(!require(rhdf5)){
    print("please install rhdf5")}
  else{
    #library(rhdf5)
    #read the interactions matrix
    #file<-"/lustre/user/liclab/licontinent/wxm/wxm/hic1/GSE105465_ENCFF762DPV_chromatin_interactions_hg19.h5"
    interactions<-h5read(file,"interactions")
    #print(interactions[1:2,1:4])
    #deal with the NaN is in plotTAD
    #read the bin positions and get the resolution
    bins<-h5read(file,"bin_positions")
    resolution<-bins[2,2]-bins[2,1]
    print(paste("this hic data's resolution is",resolution/1000,"kb"))
    #read the chrs
    chrs<-h5read(file,"chrs")
    #read the  chr_bins_range<-h5read
    chr_bins_range<-h5read(file,"chr_bin_range")
    from1<-from/resolution
    to1<-to/resolution
    if(class(chr)=="numeric"){
      from1<-from1+chr_bins_range[1,chr]
      to1<-to1+chr_bins_range[1,chr]
      num<-chr
    }
    else{
      for(i in 1:length(chrs)){
        if(chr==chrs[i]){
          from1<-from1+chr_bins_range[1,i]
          to1<-to1+chr_bins_range[1,i]
          num<-i
        }
      }
    }
    if((to1>chr_bins_range[2,num])){
      print("error: out of range of chromosome",chr)
      print("please check the range of the chromosome")
      colnames(chr_bins_range)<-chrs
      print(paste("----------------------------------"))
      print("the chr bin range and chr range:")
      print("the chr bin range is ")
      print(chr_bins_range)
      chr_bins_range<-as.numeric( chr_bins_range)
      print("the chr range is ")
      print(chr_bins_range*resolution)
      print(paste("----------------------------------"))
      print("you can also use the function see_hic_h5() to get more information")
      print(paste("----------------------------------"))
      return(NULL)
    }
    plothicdata<-interactions[from1:to1,from1:to1]
    # 参数说明
    ## input.matrix: 输入矩阵
    ## col.min="red": 热图颜色，支持字符或者是RGB参数za
    ## col.max="red": 热图颜色，支持字符或者是RGB参数
    ## col.boundary=0: 两种颜色的分界线，也即白色对应的值
    ## ylim=NULL: 热图的y轴范围，1个bin的长度是1，方便截短
    ## minBound=0: 使用分位数来修正最小值
    ## maxBound=0.95: 使用分位数修正最大值，默认matrix中的95%分位数为最大值
    ## mat.upper=T: 如果为TRUE画上半部分，为FALSE画下半部分; 只有在mat.part=T时有效
    ## mat.part=T: 如果为TRUE则画三角，如果为FLASE则画倒置的matrix
    input.matrix = as.matrix(plothicdata)
    for(i in 1:dim(input.matrix)[1]){
      for(j in 1:dim(input.matrix)[1]){
        if(is.nan(input.matrix[i,j])){
          input.matrix[i,j]=0
        }
        if(is.na(input.matrix[i,j]))
          input.matrix[i,j]=0
      }
    }
    mat.nrow = dim(input.matrix)[1]#dim会产生行和列的参数，第一个为几行
    
    ## color type 1：single color
    ## color type 2：blue white red black
    
    ## 生成矩阵需要的颜色
    #按照列取出对角线下方（不含对角线）的值，并转化为向量
    mat.lower.tri = as.vector(input.matrix[lower.tri(input.matrix,diag = F)])
    #取出上下分位数
    mat.quantile = quantile(as.vector(input.matrix[lower.tri(input.matrix,diag = T)]),prob=c(minBound,maxBound))
    
    if(col.boundary>mat.quantile[2] | col.boundary<mat.quantile[1]){
      print("Error! col.boundary have to smaller than matrix maxBound value!")
      return(NULL)
    }
    
    ## fix value,too large or too small
    #此时他修改了对角线和对角线下的数据
    mat.lower.tri[mat.lower.tri < mat.quantile[1]] = mat.quantile[1]
    mat.lower.tri[mat.lower.tri > mat.quantile[2]] = mat.quantile[2]
    mat.diagnal.value = diag(input.matrix)
    mat.diagnal.value[mat.diagnal.value < mat.quantile[1]] = mat.quantile[1]
    mat.diagnal.value[mat.diagnal.value > mat.quantile[2]] = mat.quantile[2]
    
    ## create color vector
    #初始化向量，各个值都为0
    mat.lower.tri.col_alpha = rep(0,length(mat.lower.tri))
    mat.diagnal.col_alpha = rep(0,length(mat.diagnal.value))
    #ceiling 指小数进1
    #255是线性变化的加权值
    mat.lower.tri.col_alpha[mat.lower.tri>=col.boundary] = ceiling((mat.lower.tri[mat.lower.tri>=col.boundary] - col.boundary) / (mat.quantile[2] - col.boundary) * 255)
    mat.lower.tri.col_alpha[mat.lower.tri<col.boundary] = ceiling((col.boundary - mat.lower.tri[mat.lower.tri<col.boundary]) / (col.boundary - mat.quantile[1]) * 255)
    mat.diagnal.col_alpha[mat.diagnal.value>=col.boundary] = ceiling((mat.diagnal.value[mat.diagnal.value>=col.boundary] - col.boundary) / (mat.quantile[2] - col.boundary) * 255)
    mat.diagnal.col_alpha[mat.diagnal.value<col.boundary] = ceiling((col.boundary - mat.diagnal.value[mat.diagnal.value<col.boundary]) / (col.boundary - mat.quantile[1]) * 255)
    #初始化颜色矩阵，默认##FFFFFF为白色，相当于背景板
    mat.lower.tri.color = rep("##FFFFFF",length(mat.lower.tri))
    mat.diagnal.color = rep("##FFFFFF",length(mat.diagnal.value))
    #rgb()生成颜色编码
    #col2rgb()将颜色转为rgb色值。alpha:设定透明度
    mat.lower.tri.color[mat.lower.tri>=col.boundary] = rgb(t(col2rgb(col.max)),alpha = mat.lower.tri.col_alpha[mat.lower.tri>=col.boundary],maxColorValue = 255)
    mat.lower.tri.color[mat.lower.tri<col.boundary] = rgb(t(col2rgb(col.min)),alpha = mat.lower.tri.col_alpha[mat.lower.tri<col.boundary],maxColorValue = 255)
    mat.diagnal.color[mat.diagnal.value>=col.boundary] = rgb(t(col2rgb(col.max)),alpha = mat.diagnal.col_alpha[mat.diagnal.value>=col.boundary],maxColorValue = 255)
    mat.diagnal.color[mat.diagnal.value<col.boundary] = rgb(t(col2rgb(col.min)),alpha = mat.diagnal.col_alpha[mat.diagnal.value<col.boundary],maxColorValue = 255)
    
    
    
    # 根据矩阵的行列，生成x1,y1的坐标 全用的行数
    x1.part = c(1:(mat.nrow -1))
    y1.part = c(1:(mat.nrow -1))
    x1 = x1.part
    y1 = y1.part
    ##  ##这个向量的变化很不错，可以以三维矩阵为例搞清楚它
    for(i in c(1:(mat.nrow -1))){
      x1.part = x1.part[-length(x1.part)] + 2
      y1.part = y1.part[-length(y1.part)]
      
      x1 = c(x1,x1.part)
      y1 = c(y1,y1.part)
    }
    
    x1# 对角线的x1
    x1.diagonal = c(0:(mat.nrow-1)) * 2
    y1.diagonal = rep(0,mat.nrow)
    
    if(mat.part==T){#画三角
      # plot matrix upper OR lower only.
      x1 = c(x1,x1.diagonal)
      y1 = c(y1,y1.diagonal)  
      # 颜色向量也不同
      mat.color = c(mat.lower.tri.color,mat.diagnal.color)
    }
    else if(mat.part==F){#画矩阵
      x1 = c(x1,x1.diagonal,x1)
      y1 = c(y1,y1.diagonal,-y1)  
      mat.color = c(mat.lower.tri.color,mat.diagnal.color,mat.lower.tri.color)
    }
    
    
    # 根据x1 y1 生成剩余坐标         #x2,y2
    x2 = x1 + 1                    
    x3 = x2 + 1                 #x1,y1    #x3,y3
    x4 = x2
    y2 = y1 + 1                      #x4,y4
    y3 = y1
    y4 = y1 - 1
    # 需要按照x1,x2,x3,x4,NA的格式进行生成x vector否则会把所有的点连在一起
    NA_vector = rep(NA,length(x1))
    # 生成按照x11,x12,x13,x14,NA,x21,x22,x23,x24,NA...排列的x与y
    x_matrix = matrix(c(x1,x2,x3,x4,NA_vector),ncol = 5)
    y_matrix = matrix(c(y1,y2,y3,y4,NA_vector),ncol = 5)
    x = as.vector(t(x_matrix))
    y = as.vector(t(y_matrix))
    
    # 绘图部分
    ## 如果mat.upper 为F 则画倒置的图像
    ## 设置画布大小
    x.point = c(0,mat.nrow*2)
    if(mat.part){
      if(! mat.upper){ y = -y ; y.point = c(-mat.nrow,1)}
      else{ y.point = c(-1,mat.nrow) }
    }
    else{ 
      y.point = c(-mat.nrow,mat.nrow)
    }
    ## 设置ylim取值范围
    if(is.null(ylim)){ ylim = y.point }
    
    ## 生成画布  ####xaxt="n",yaxt="n"禁用xy刻度线
    plot(x.point,y.point,ylim=ylim,type="n",frame.plot =F,xaxt="n",yaxt="n",xlab="",ylab="")
    ## 画三角矩阵
    polygon(x,y,col = mat.color,border = F)
    if(title){title(paste("plot hic of",chr,"from:",from,"to",to,"\nresolution:",resolution/1000,"kb"))}
  }
}




#plot inter-chromosome interactions 
ploth5_inter_hic<-function(file,chr_a,chr_b,col.min = "red",col.max = "red",col.boundary = 0,ylim=NULL,minBound=0,maxBound=0.95,mat.upper=T,mat.part=T){
  
  
  #file :h5 file from NCBI format is above
  #chr plot which chromosome
  #from to position such as chr1,1000000,1500000
  if(!require(rhdf5)){
    return("please install rhdf5")}
  else{
    library(rhdf5)
    #read the interactions matrix
    #file<-"/lustre/user/liclab/licontinent/wxm/wxm/hic1/GSE105465_ENCFF762DPV_chromatin_interactions_hg19.h5"
    interactions<-h5read(file,"interactions")
    #deal with the NaN is in plotTAD
    #read the bin positions and get the resolution
    bins<-h5read(file,"bin_positions")
    resolution<-bins[2,2]-bins[2,1]
    print(paste("this hic data's resolution is",resolution/1000,"kb"))
    #read the chrs
    chrs<-h5read(file,"chrs")
    #read the  chr_bins_range<-h5read
    #chr_bins_range<-h5read(file,"chr_bin_range")
    #from<-from/resolution
    #to<-to/resolution
    if(class(chr_a)=="numeric"){
      chra_begin<-chr_bins_range[1,chr_a]
      chra_end<-chr_bins_range[2,chr_a]
      chrb_begin<-chr_bins_range[1,chr_b]
      chrb_end<-chr_bins_range[2,chr_b]
    }
    else{
      for(i in 1:length(chrs)){
        if(chr_a==chrs[i]){
          chra_begin<-chr_bins_range[1,i]
          chra_end<-chr_bins_range[2,i]
        }
        if(chr_b==chrs[i]){
          chrb_begin<-chr_bins_range[1,i]
          chrb_end<-chr_bins_range[2,i]
        }
      }
    }
    mat_aa<-interactions[chra_begin:chra_end,chra_begin:chra_end]
    mat_ab<-interactions[chra_begin:chra_end,chrb_begin:chrb_end]
    mat_ba<-interactions[chrb_begin:chrb_end,chra_begin:chra_end]
    mat_bb<-interactions[chrb_begin:chrb_end,chrb_begin:chrb_end]
    plothicdata<-rbind(cbind(mat_aa,mat_ab),cbind(mat_ba,mat_bb))
    # 参数说明
    ## input.matrix: 输入矩阵
    ## col.min="red": 热图颜色，支持字符或者是RGB参数za
    ## col.max="red": 热图颜色，支持字符或者是RGB参数
    ## col.boundary=0: 两种颜色的分界线，也即白色对应的值
    ## ylim=NULL: 热图的y轴范围，1个bin的长度是1，方便截短
    ## minBound=0: 使用分位数来修正最小值
    ## maxBound=0.95: 使用分位数修正最大值，默认matrix中的95%分位数为最大值
    ## mat.upper=T: 如果为TRUE画上半部分，为FALSE画下半部分; 只有在mat.part=T时有效
    ## mat.part=T: 如果为TRUE则画三角，如果为FLASE则画倒置的matrix
    input.matrix = as.matrix(plothicdata)
    for(i in 1:dim(input.matrix)[1]){
      for(j in 1:dim(input.matrix)[1]){
        if(is.nan(input.matrix[i,j])){
          input.matrix[i,j]=0
        }
        if(is.na(input.matrix[i,j]))
          input.matrix[i,j]=0
      }
    }
    mat.nrow = dim(input.matrix)[1]#dim会产生行和列的参数，第一个为几行
    
    ## color type 1：single color
    ## color type 2：blue white red black
    
    ## 生成矩阵需要的颜色
    #按照列取出对角线下方（不含对角线）的值，并转化为向量
    mat.lower.tri = as.vector(input.matrix[lower.tri(input.matrix,diag = F)])
    #取出上下分位数
    mat.quantile = quantile(as.vector(input.matrix[lower.tri(input.matrix,diag = T)]),prob=c(minBound,maxBound))
    
    if(col.boundary>mat.quantile[2] | col.boundary<mat.quantile[1]){
      print("Error! col.boundary have to smaller than matrix maxBound value!")
      return(NULL)
    }
    
    ## fix value,too large or too small
    #此时他修改了对角线和对角线下的数据
    mat.lower.tri[mat.lower.tri < mat.quantile[1]] = mat.quantile[1]
    mat.lower.tri[mat.lower.tri > mat.quantile[2]] = mat.quantile[2]
    mat.diagnal.value = diag(input.matrix)
    mat.diagnal.value[mat.diagnal.value < mat.quantile[1]] = mat.quantile[1]
    mat.diagnal.value[mat.diagnal.value > mat.quantile[2]] = mat.quantile[2]
    
    ## create color vector
    #初始化向量，各个值都为0
    mat.lower.tri.col_alpha = rep(0,length(mat.lower.tri))
    mat.diagnal.col_alpha = rep(0,length(mat.diagnal.value))
    #ceiling 指小数进1
    #255是线性变化的加权值
    mat.lower.tri.col_alpha[mat.lower.tri>=col.boundary] = ceiling((mat.lower.tri[mat.lower.tri>=col.boundary] - col.boundary) / (mat.quantile[2] - col.boundary) * 255)
    mat.lower.tri.col_alpha[mat.lower.tri<col.boundary] = ceiling((col.boundary - mat.lower.tri[mat.lower.tri<col.boundary]) / (col.boundary - mat.quantile[1]) * 255)
    mat.diagnal.col_alpha[mat.diagnal.value>=col.boundary] = ceiling((mat.diagnal.value[mat.diagnal.value>=col.boundary] - col.boundary) / (mat.quantile[2] - col.boundary) * 255)
    mat.diagnal.col_alpha[mat.diagnal.value<col.boundary] = ceiling((col.boundary - mat.diagnal.value[mat.diagnal.value<col.boundary]) / (col.boundary - mat.quantile[1]) * 255)
    #初始化颜色矩阵，默认##FFFFFF为白色，相当于背景板
    mat.lower.tri.color = rep("##FFFFFF",length(mat.lower.tri))
    mat.diagnal.color = rep("##FFFFFF",length(mat.diagnal.value))
    #rgb()生成颜色编码
    #col2rgb()将颜色转为rgb色值。alpha:设定透明度
    mat.lower.tri.color[mat.lower.tri>=col.boundary] = rgb(t(col2rgb(col.max)),alpha = mat.lower.tri.col_alpha[mat.lower.tri>=col.boundary],maxColorValue = 255)
    mat.lower.tri.color[mat.lower.tri<col.boundary] = rgb(t(col2rgb(col.min)),alpha = mat.lower.tri.col_alpha[mat.lower.tri<col.boundary],maxColorValue = 255)
    mat.diagnal.color[mat.diagnal.value>=col.boundary] = rgb(t(col2rgb(col.max)),alpha = mat.diagnal.col_alpha[mat.diagnal.value>=col.boundary],maxColorValue = 255)
    mat.diagnal.color[mat.diagnal.value<col.boundary] = rgb(t(col2rgb(col.min)),alpha = mat.diagnal.col_alpha[mat.diagnal.value<col.boundary],maxColorValue = 255)
    
    # 根据矩阵的行列，生成x1,y1的坐标 全用的行数
    x1.part = c(1:(mat.nrow -1))
    y1.part = c(1:(mat.nrow -1))
    x1 = x1.part
    y1 = y1.part
    ##  ##这个向量的变化很不错，可以以三维矩阵为例搞清楚它
    for(i in c(1:(mat.nrow -1))){
      x1.part = x1.part[-length(x1.part)] + 2
      y1.part = y1.part[-length(y1.part)]
      
      x1 = c(x1,x1.part)
      y1 = c(y1,y1.part)
    }
    
    x1# 对角线的x1
    x1.diagonal = c(0:(mat.nrow-1)) * 2
    y1.diagonal = rep(0,mat.nrow)
    
    if(mat.part==T){#画三角
      # plot matrix upper OR lower only.
      x1 = c(x1,x1.diagonal)
      y1 = c(y1,y1.diagonal)  
      # 颜色向量也不同
      mat.color = c(mat.lower.tri.color,mat.diagnal.color)
    }
    else if(mat.part==F){#画矩阵
      x1 = c(x1,x1.diagonal,x1)
      y1 = c(y1,y1.diagonal,-y1)  
      mat.color = c(mat.lower.tri.color,mat.diagnal.color,mat.lower.tri.color)
    }
    
    
    # 根据x1 y1 生成剩余坐标         #x2,y2
    x2 = x1 + 1                    
    x3 = x2 + 1                 #x1,y1    #x3,y3
    x4 = x2
    y2 = y1 + 1                      #x4,y4
    y3 = y1
    y4 = y1 - 1
    # 需要按照x1,x2,x3,x4,NA的格式进行生成x vector否则会把所有的点连在一起
    NA_vector = rep(NA,length(x1))
    # 生成按照x11,x12,x13,x14,NA,x21,x22,x23,x24,NA...排列的x与y
    x_matrix = matrix(c(x1,x2,x3,x4,NA_vector),ncol = 5)
    y_matrix = matrix(c(y1,y2,y3,y4,NA_vector),ncol = 5)
    x = as.vector(t(x_matrix))
    y = as.vector(t(y_matrix))
    
    # 绘图部分
    ## 如果mat.upper 为F 则画倒置的图像
    ## 设置画布大小
    x.point = c(0,mat.nrow*2)
    if(mat.part){
      if(! mat.upper){ y = -y ; y.point = c(-mat.nrow,1)}
      else{ y.point = c(-1,mat.nrow) }
    }
    else{ 
      y.point = c(-mat.nrow,mat.nrow)
    }
    ## 设置ylim取值范围
    if(is.null(ylim)){ ylim = y.point }
    
    ## 生成画布  ####xaxt="n",yaxt="n"禁用xy刻度线
    plot(x.point,y.point,ylim=ylim,type="n",frame.plot =F,xaxt="n",yaxt="n",xlab="",ylab="")
    ## 画三角矩阵
    polygon(x,y,col = mat.color,border = F)
    title(paste("the inter-chro hic plot of",chr_a,"and",chr_b))
  }
}
