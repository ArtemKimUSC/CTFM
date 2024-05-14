library(susieR)
args = commandArgs(trailingOnly=TRUE)

wd=args[1]
trait=args[2]


### read correlation matrix EUR ###
cor.res=read.table(paste0(wd,'/data/SLDSC_default/R.txt'), sep='\t', header=T) # read the correlation matrix
cor.res=as.matrix(cor.res)
annots=read.table(paste0(wd,'/data/SLDSC_default/annots.txt'), header=T, sep='\t')
cor.res=cor.res[rownames(cor.res) %in% annots$Name, colnames(cor.res) %in% annots$Name]


all_sldsc=data.frame()
all_pips=data.frame()

i=trait

    data=data.frame()

    for (j in 0:5) {
    tmp=read.table(paste0(wd,'/out/SLDSC/',i,'/',i,'.',j,'.cell_type_results.txt'), sep='\t', header=T)
    data=rbind(data,tmp)
    }
    rm(tmp)




#######
data=data[!duplicated(data$Name),]

data$Name=gsub('ABC_old_','',data$Name)
data$Name=gsub('-','.',data$Name)
data$Name=gsub('Zhang_Zhang_','Zhang_',data$Name)
data$Name=gsub('ENCODE_','',data$Name)

if (all(colnames(cor.res) %in% data$Name)==TRUE) {print('match OK')} else {print('ERROR match'); stop()}

data=data[data$Name %in% rownames(cor.res),]
if (nrow(data)==927) { print('OK');print(i) } else {print('ERROR: Number of annots is not 927, stopping...');next}


data=data[match(rownames(cor.res), data$Name), ]


data$Zscore=data$Coefficient/data$Coefficient_std_error
data$TRAIT=i



#out1 => for pvalues SLDSC
all_sldsc=rbind(all_sldsc,data)


#remove negZ for susie
data=data[data$Zscore>0,]
cor.res_temp=cor.res[rownames(cor.res) %in% data$Name, colnames(cor.res) %in% data$Name]
data=data[match(rownames(cor.res_temp), data$Name), ]

#Run SuSiE
rss=susie_rss(R=cor.res_temp, z=data$Zscore, estimate_residual_variance = F, L=10, n=1190321)

#out2 => pips
pip=as.data.frame(rss$pip)
pip$Name=data$Name
pip$Zscore=data$Zscore
pip$Coeff=data$Coefficient
pip$TRAIT=i
pip$INDEX=1:nrow(pip)
colnames(pip)[1]='PIP'
all_pips=rbind(all_pips,pip)


#out3 => cs
all_cs=data.frame()

a=susie_get_cs(rss,Xcorr=cor.res_temp,coverage=0.96,min_abs_corr = 0.3) # maybe use Xcorr to filter out based on correlation threshold
                    #Xcorr : p by p matrix of correlations between variables (covariates). When provided, it will be used to remove CSs whose minimum correlation among variables is smaller than min_abs_corr.

if (length(a$cs)==0) {print(paste0(i,' no CS with SuSiE purity threshold')); 


} else for (j in 1:length(a$cs)) {
    cs=data.frame( 'CS_Num'=names(a$cs)[j],'Size'=length(a$cs[[j]]),'annots'=a$cs[[j]], 'TRAIT'=trait)
    cs$annot_name=pip$Name[match(cs$annots,pip$INDEX)]
    cs$INDEX=pip$INDEX[match(cs$annot_name,pip$Name)]
    cs$PIP=pip$PIP[match(cs$annots,pip$INDEX)]
    cs$Zscore=pip$Zscore[match(cs$annots,pip$INDEX)]
    cs$Alpha=rss$alpha[j,cs$INDEX]
    cs$Type='SuSiE'
    all_cs=rbind(all_cs,cs)
}

print(paste0('Number of SuSiE CS: ', nrow(table(all_cs$CS_Num))))

if (nrow(all_cs) >=1) {all_cs$S=0}


L=data.frame()
cs_alpha=data.frame()
S_all=data.frame()

if (nrow(rss$alpha)>=10) {

for (i in 1:10) { S=sum(rss$alpha[i,] * log(rss$alpha[i,] * length(rss$alpha[i,]))); if (S>3) {
l=data.frame(L=paste0('L',i)); L=rbind(L,l); 
s=data.frame(L=paste0('L',i),S=S); S_all=rbind(S_all,s) } 
}
}



if (nrow(L) > 0) {

for (i in 1:nrow(L)) {if (!(L$L[i] %in% names(a$cs))) {
    b=as.data.frame(rss$alpha[i,])
    b$annot_name=row.names(b)
    b=b[order(b[,1],decreasing=T),]
    SUM=0
    j=1
    cs=data.frame()
    while (SUM < 0.9) {
        cs=rbind(cs,b[j,])
        SUM=SUM+b[j,1]
        j=j+1

    }
    cs$CS_Num=L$L[i]
    cs$Size=nrow(cs)
    cs$PIP=pip$PIP[match(cs$annot_name,pip$Name)]
    cs$TRAIT=trait
    cs$Zscore=pip$Zscore[match(cs$annot_name,pip$Name)]
    cs$Type='Alpha'
    colnames(cs)[1]='Alpha'
    cs$S=S_all$S[match(cs$CS_Num,S_all$L)]
    cs$INDEX=pip$INDEX[match(cs$annot_name,pip$Name)]
    cs$annots=pip$INDEX[match(cs$annot_name,pip$Name)]
    if (nrow(all_cs > 0)) { cs=cs[,colnames(all_cs)] }
 
    cs_alpha=rbind(cs_alpha,cs)

    } else { print(paste(L$L[i],' already in Susie CS')) }  }
    } else {  
        print('no alpha CS')
}

print(paste0('Number of alpha CS: ', nrow(table(cs_alpha$CS_Num))))

all_cs=rbind(all_cs,cs_alpha)




print(paste0("TRAIT:",trait))
print(paste0("number of total CS:",nrow(table(all_cs$CS_Num))))


write.table(all_pips, file=paste0(wd,'/out/susie/',trait,'_CTFM_pips.txt'),sep='\t', quote=F,row.names=F)
write.table(all_sldsc, file=paste0(wd,'/out/susie/',trait,'_CTFM_sldsc.txt'),sep='\t', quote=F,row.names=F)
write.table(all_cs, file=paste0(wd,'/out/susie/',trait,'_CTFM_CS.txt'),sep='\t', quote=F,row.names=F)


print(paste0('DONE ',trait))







