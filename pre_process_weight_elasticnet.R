
setwd("/home/panwei/wuxx0845/TWAS_casual/WEIGHTS/GTEx.ALL")

library(data.table)
tissue = list.dirs('.', recursive=FALSE)

tissue = gsub("./", "", tissue)


final.dat = NULL
for( i in 1:length(tissue)) {
    if(i%in%c(1,45,46,47)) {
        tmp.pos = read.table(paste(tissue[i],"/",tissue[i],".pos",sep =""),header = TRUE)
        tmp.pos[,1] = as.character(tmp.pos[,1])
        tmp.pos$id2 = paste(gsub("\\..*","",tmp.pos[,1]),tmp.pos[,2],sep=".")
        rownames(tmp.pos) = tmp.pos$id2
        
        tmp.profile = read.table(paste(tissue[i],"/",tissue[i],".profile",sep =""),header = TRUE)
        rownames(tmp.profile) = tmp.profile$id
        
        tmp.profile = tmp.profile[tmp.pos$id2,]
        tmp = cbind(tmp.pos,tmp.profile)
        final.dat = rbind(final.dat, tmp)
    } else {
        tmp.pos = read.table(paste(tissue[i],"/",tissue[i],".pos",sep =""),header = TRUE)
        tmp.pos[,1] = as.character(tmp.pos[,1])
        
        tmp.name = gsub("^.*/","",tmp.pos[,1])
        tmp.name = gsub(".wgt.RDat","",tmp.name)
        tmp.pos$id2 = tmp.name
        rownames(tmp.pos) = tmp.pos$id2
        
        tmp.profile = read.table(paste(tissue[i],"/",tissue[i],".profile",sep =""),header = TRUE)
        rownames(tmp.profile) = tmp.profile$id
        
        tmp.profile = tmp.profile[tmp.pos$id2,]
        tmp = cbind(tmp.pos,tmp.profile)
        

        #colSums(is.na(tmp))
        #dim(tmp)
        final.dat = rbind(final.dat, tmp)
    }
    
}

final.dat2 = final.dat[,c(-6)]

gene.list = unique(final.dat2[,"ID"])

info.save = NULL

length(gene.list)

time.start = proc.time()
for(i in 1:length(gene.list)) { #length(gene.list)
    tmp = final.dat2[final.dat2[,"ID"]==gene.list[i],]
    
    #dim(tmp)
    rownames(tmp) = 1:dim(tmp)[1]
    
    # tmp$max.r2 <- apply(tmp[,12:14], 1, function(x) max(x,na.rm=TRUE))
    # tmp$max.r2.2 <- apply(tmp[,11:14], 1, function(x) max(x,na.rm=TRUE))
    
    tmp$enet.r2[!is.finite(tmp$enet.r2)] = -1
    
    tmp = tmp[tmp$enet.r2==max(tmp$enet.r2),]
    
    if (tmp$enet.r2 == -1) {
        next
    }
    
    tmp = tmp[1,]
    
    sub.dir = strsplit(as.character(tmp[1]),"/")[[1]]
    
    tmp$WGT2 = paste("twas.com.net.weight/",sub.dir[2],sep="")
    
    from= paste(sub.dir[1],"/",tmp[1],sep="")
    to = paste("/home/panwei/wuxx0845/TWAS_casual/WEIGHTS/twas.com.net.weight/",sub.dir[2],sep="")
    file.copy(from,to)
    info.save = rbind(info.save,tmp)
}

saveRDS(info.save,"/home/panwei/wuxx0845/TWAS_casual/WEIGHTS/twas.com.net.info.save.rds")

info.pos = info.save[,c("WGT2","ID","CHR","P0","P1")]
colnames(info.pos) = c("WGT","ID","CHR","P0","P1")

info.pos = info.pos[order(info.pos[,3]),]
write.table(info.pos,"/home/panwei/wuxx0845/TWAS_casual/WEIGHTS/twas.com.net.weight.pos",quote = FALSE)

info.profile = info.save[,c("id","nsnps","hsq","hsq.se","hsq.pv","top1.r2","blup.r2","enet.r2","bslmm.r2","lasso.r2","top1.pv","blup.pv","enet.pv","bslmm.pv","lasso.pv")]

write.table(info.profile,"/home/panwei/wuxx0845/TWAS_casual/WEIGHTS/twas.com.net.weight.profile",quote = FALSE)


# generate METSIM based weights
info.save = readRDS("/home/panwei/wuxx0845/TWAS_casual/WEIGHTS/twas.com.net.info.save.rds")


info.pos = info.save[,c("WGT2","ID","CHR","P0","P1")]
colnames(info.pos) = c("WGT","ID","CHR","P0","P1")

rownames(info.pos) = info.pos[,2]

wgtlist = read.table("/home/panwei/wuxx0845/TWAS_casual/WEIGHTS/METSIM.ADIPOSE.RNASEQ.pos",head=T,as.is=T)

rownames(wgtlist) = wgtlist[,2]

info.pos2 = info.pos[!info.pos[,2] %in%wgtlist[,2],]
info.pos3 = rbind(info.pos2,wgtlist)

info.pos3 = info.pos3[order(info.pos3[,3]),]
write.table(info.pos3,"/home/panwei/wuxx0845/TWAS_casual/WEIGHTS/METSIM.twas.com.net.weight.pos",quote = FALSE)


# generate Brain based weights
info.save = readRDS("/home/panwei/wuxx0845/TWAS_casual/WEIGHTS/twas.com.net.info.save.rds")


info.pos = info.save[,c("WGT2","ID","CHR","P0","P1")]
colnames(info.pos) = c("WGT","ID","CHR","P0","P1")

rownames(info.pos) = info.pos[,2]

wgtlist = read.table("/home/panwei/wuxx0845/TWAS_casual/WEIGHTS/CMC.BRAIN.RNASEQ.pos",head=T,as.is=T)

rownames(wgtlist) = wgtlist[,2]

info.pos2 = info.pos[!info.pos[,2] %in%wgtlist[,2],]
info.pos3 = rbind(info.pos2,wgtlist)

info.pos3 = info.pos3[order(info.pos3[,3]),]

write.table(info.pos3,"/home/panwei/wuxx0845/TWAS_casual/WEIGHTS/CMC.twas.com.net.weight.pos",quote = FALSE)
