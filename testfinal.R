data = read.csv("Network10.txt",
                sep=",",header = T)
update_original = read.csv("matlab.csv",
                  header=T,sep=",")
names(update_original) = c(data[,1])
update_original = as.matrix(update_original)



binary = data.frame(
  c(      rep(0,512),rep(1,512)),
  c(rep(c(rep(0,256),rep(1,256)),2)),
  c(rep(c(rep(0,128),rep(1,128)),4)),
  c(rep(c(rep(0, 64),rep(1, 64)),4)),
  c(rep(c(rep(0, 32),rep(1, 32)),4)),
  c(rep(c(rep(0, 16),rep(1, 16)),4)),
  c(rep(c(rep(0,  8),rep(1,  8)),4)),
  c(rep(c(rep(0,  4),rep(1,  4)),4)),
  c(rep(c(rep(0,  2),rep(1,  2)),4)),
  c(rep(c(rep(0,  1),rep(1,  1)),4)))

names(binary) = c(data[,1])
binary = as.matrix(binary)

###############################
range=30
scheme=10
n0=1
n1=1000
mat = array(c(0,0,0),dim = c(1,3))
k=1

nrow(update)

for(i in seq(1,9,1000)){
  n0=i
  n1=i+8
  
  update = update_original[which(apply(update_original,1,
                                       function(x) max(x[1:10])==scheme)==TRUE),]
  final_bin = array(dim = c(1024,10,length(n0:n1),range+1))
  final_bin[,,,1] <- binary
  count = 1
  attr1 = c(0,0,0,0,0,0,0,0,0,0)
  attr2 = c(1,0,1,0,0,0,0,0,1,1)
  attr3 = c(1,1,0,1,1,1,1,1,0,0)
  

  
  for(j in 1:range){
    
    #pre-loading conditions
    condition1.1 = update[n0:n1,1] == count
    condition1.2 = final_bin[,2,,j] == 1 |
      final_bin[,4,,j] == 1 |
      final_bin[,9,,j] == 1
    
    condition2.1 = update[n0:n1,2] == count
    condition2.2 = final_bin[,1,,j] == 1 &
      (final_bin[,3,,j] == 0 & final_bin[,9,,j] == 0)
    
    condition3.1 = update[n0:n1,3] == count
    condition3.2 = (final_bin[,1,,j] == 1 |
                      (final_bin[,9,,j] == 1 & final_bin[,10,,j] == 1)) &
      (final_bin[,2,,j] == 0 & final_bin[,5,,j] == 0)
    
    condition4.1 = update[n0:n1,4] == count
    condition4.2 = (final_bin[,4,,j] == 1 |
                      final_bin[,5,,j] == 1) &
      final_bin[,10,,j] == 0
    
    condition5.1 = update[n0:n1,5] == count
    condition5.2 = ((final_bin[,5,,j] == 1 & final_bin[,6,,j] == 1) |
                      (final_bin[,1,,j] == 1 & final_bin[,2,,j] == 1 & final_bin[,4,,j] == 1) |
                      final_bin[,4,,j] == 1 |
                      final_bin[,7,,j] == 1) &
      (final_bin[,9,,j] == 0 | final_bin[,10,,j] == 0)
    
    condition6.1 = update[n0:n1,6] == count
    condition6.2 = final_bin[,5,,j] == 1 &
      final_bin[,9,,j] == 0
    
    condition7.1 = update[n0:n1,7] == count
    condition7.2 = final_bin[,5,,j] == 1
    
    condition8.1 = update[n0:n1,8] == count
    condition8.2 = final_bin[,4,,j] == 1 |
      final_bin[,5,,j] == 1
    
    condition9.1 = update[n0:n1,9] == count
    condition9.2 = (final_bin[,1,,j] == 1 |
                      final_bin[,3,,j] == 1 |
                      final_bin[,9,,j] == 1 |
                      final_bin[,10,,j] == 1) &
      (final_bin[,6,,j] == 0 & final_bin[,8,,j] == 0)
    
    condition10.1 = update[n0:n1,10] == count
    condition10.2 = (final_bin[,9,,j] == 1 | final_bin[,10,,j] == 1) &
      (final_bin[,4,,j] == 0 |
         (final_bin[,1,,j] == 0 & final_bin[,2,,j] == 0 & final_bin[,5,,j] == 0) |
         final_bin[,8,,j] == 0)
    
    
    final_bin[,,,j+1] <- array(mapply(updt,length(n0:n1)),dim=c(1024,10,length(n0:n1)))
    
    #updating cycle
    if(count==max(update[k,])){
      count=1
    }else{
      count=count+1
    }
    
  }
  
  mat[1,1] <- sum(apply(final_bin[,,,range+1],3,
                        function (x) nrow(subset(x,x[,1]==0 &
                                                   x[,2]==0 &
                                                   x[,3]==0 &
                                                   x[,4]==0 &
                                                   x[,5]==0 &
                                                   x[,6]==0 &
                                                   x[,7]==0 &
                                                   x[,8]==0 &
                                                   x[,9]==0 &
                                                   x[,10]==0))))+mat[1,1]
                                                                            
  mat[1,2] <- sum(apply(final_bin[,,,range+1],3,
                        function (x) nrow(subset(x,x[,1]==1 &
                                                   x[,2]==0 &
                                                   x[,3]==1 &
                                                   x[,4]==0 &
                                                   x[,5]==0 &
                                                   x[,6]==0 &
                                                   x[,7]==0 &
                                                   x[,8]==0 &
                                                   x[,9]==1 &
                                                   x[,10]==1))))+mat[1,2]
                                                                            
  mat[1,3] <- sum(apply(final_bin[,,,range+1],3,
                        function (x) nrow(subset(x,x[,1]==1 &
                                                   x[,2]==1 &
                                                   x[,3]==0 &
                                                   x[,4]==1 &
                                                   x[,5]==1 &
                                                   x[,6]==1 &
                                                   x[,7]==1 &
                                                   x[,8]==1 &
                                                   x[,9]==0 &
                                                   x[,10]==0))))+mat[1,3]
                                                                            
  print(mat)
    
}


{
  node1 = function(x,y,z){
    
    if(condition1.1[x]){
      if(condition1.2[y]){
        1
      }else{
        0
      }
    }else{
      final_bin[y,1,x,z]
      
    }
    
  }
  
  node2 = function(x,y,z){
    
    if(condition2.1[x]){
      if(condition2.2[y]){
        1
      }else{
        0
      }
    }else{
      final_bin[y,2,x,z]
    }
    
  }
  
  node3 = function(x,y,z){
    
    if(condition3.1[x]){
      if(condition3.2[y]){
        1
      }else{
        0
      }
    }else{
      final_bin[y,3,x,z]
    }
    
  }
  
  node4 = function(x,y,z){
    
    if(condition4.1[x]){
      if(condition4.2[y]){
        1
      }else{
        0
      }
    }else{
      final_bin[y,4,x,z]
    }
    
  }
  
  node5 = function(x,y,z){
    
    if(condition5.1[x]){
      if(condition5.2[y]){
        1
      }else{
        0
      }
    }else{
      final_bin[y,5,x,z]
    }
    
  }
  
  node6 = function(x,y,z){
    
    if(condition6.1[x]){
      if(condition6.2[y]){
        1
      }else{
        0
      }
    }else{
      final_bin[y,6,x,z]
    }
    
  }
  
  node7 = function(x,y,z){
    
    if(condition7.1[x]){
      if(condition7.2[y]){
        1
      }else{
        0
      }
    }else{
      final_bin[y,7,x,z]
    }
    
  }
  
  node8 = function(x,y,z){
    
    if(condition8.1[x]){
      if(condition8.2[y]){
        1
      }else{
        0
      }
    }else{
      final_bin[y,8,x,z]
    }
    
  }
  
  node9 = function(x,y,z){
    
    if(condition9.1[x]){
      if(condition9.2[y]){
        1
      }else{
        0
      }
    }else{
      final_bin[y,9,x,z]
    }
    
  }
  
  node10 = function(x,y,z){
    
    if(condition10.1[x]){
      if(condition10.2[y]){
        1
      }else{
        0
      }
    }else{
      final_bin[y,10,x,z]
    }
    
  }
  
  updt = function(k){
    temp = array(c(mapply(node1,x=k,y=1:1024,z=j),
                   mapply(node2,x=k,y=1:1024,z=j),
                   mapply(node3,x=k,y=1:1024,z=j),
                   mapply(node4,x=k,y=1:1024,z=j),
                   mapply(node5,x=k,y=1:1024,z=j),
                   mapply(node6,x=k,y=1:1024,z=j),
                   mapply(node7,x=k,y=1:1024,z=j),
                   mapply(node8,x=k,y=1:1024,z=j),
                   mapply(node9,x=k,y=1:1024,z=j),
                   mapply(node10,x=k,y=1:1024,z=j)),dim = c(1024,10))
    return(temp)
    
  }
  
}



binary_update(binary = binary,
              update = update,
              range = 12,
              n0 = 1,n1 = 990,
              scheme = 2)

a=matrix(c(1,2,3,4,3,4,1,2,3,4,3,4,1,2,3,4,3,4),nrow=6)

test=array(c(rep(a,3)),dim=c(6,3,3))

length(which(apply(update_original,1,function(x) max(x[1:10])==7)==TRUE))








