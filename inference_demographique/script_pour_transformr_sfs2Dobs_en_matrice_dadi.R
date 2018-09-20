
sfs<-read.table("/media/abdou/My_disk/Documents/M2_sepmet/Stage_M2/mil_data/dadi_mil/sfs_tout/wild_cultivated_all.obs")
x<-c()
for(i in 1:ncol(sfs)){
  for(j in 1:nrow(sfs)){
    x<-paste(x,as.character(sfs[j,i]),sep=" ")
  }
}
x

write.table(x, "/media/abdou/My_disk/Documents/M2_sepmet/Stage_M2/mil_data/dadi_mil/sfs_tout/w_c_all.fs")
