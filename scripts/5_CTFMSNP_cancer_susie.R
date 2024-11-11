args = commandArgs(trailingOnly=TRUE)
library(susieR)


#INPUT = wd + rs + trait
wd=args[1]
trait=args[2]
rs=args[3]


cor.res=read.table(paste0(wd,'/data/SLDSC_cancer/R_cancer.txt'), sep='\t', header=T) # read the correlation matrix
cor.res=as.matrix(cor.res)
#annots=read.table(paste0(wd,'/data/SLDSC_cancer/annots.txt'), header=T, sep='\t')
#cor.res=cor.res[rownames(cor.res) %in% annots$Name, colnames(cor.res) %in% annots$Name]


# Open rs.out file and rename annots

SNP=read.table(rs)
if (ncol(SNP)==8) { colnames(SNP)[8] = 'annot'} else {print("Error: the rs.out file does not contain 8 columns."); stop()}

SNP$annot=gsub('.bed','',SNP$annot)
SNP$annot=gsub("(.*/\\s*(.*$))", "\\2",SNP$annot)
SNP=SNP[!(SNP$annot %in% c('DHS_core_whole','DHS_full_whole','ABC_whole','Epimap_cell_line','Epimap_invitro','DHS_full_BRAIN.backup','Epimap_primary_cell','Epimap_tissue','Epimap_whole','GTEX_whole','Zhang_whole')),]


SNP$annot=gsub('ABC_','',SNP$annot)
SNP$annot=gsub('-','.',SNP$annot)
SNP$annot=gsub('.hg19','',SNP$annot)
SNP=SNP[SNP$annot %in% annots$Name,]

if (all(SNP$annot %in% rownames(cor.res)) != TRUE) { stop('error matching annotations between SNP file and Cor res file')}
# rownames(cor.res)=gsub('\\.','_',rownames(cor.res))
# colnames(cor.res)=gsub('\\.','_',colnames(cor.res))


### Open trait S-LDSC results and extract annotations

all_pips=data.frame()
all_cs=data.frame()
data=data.frame()

i=trait

    data=data.frame()

    for (j in 0:6) {
    tmp=read.table(paste0(wd,'/out/SLDSC_cancer/',i,'/',i,'.',j,'.cell_type_results.txt'), sep='\t', header=T)
    data=rbind(data,tmp)
    }
    rm(tmp)



data$Name=gsub('ABC_old_','',data$Name)
data$Name=gsub('Zhang_Zhang_','Zhang_',data$Name)
data$Name=gsub('ENCODE_','',data$Name)
data$Name=gsub('-','.',data$Name)

#if (all(SNP$annot %in% data$Name)==TRUE & all(SNP$annot %in% annots$Name)==TRUE & all(SNP$annot %in% rownames(cor.res))==TRUE) {print('match OK')} else {print('ERROR match'); stop()}


data=data[data$Name %in% annots$Name,]

#extract annots from S-LDSC and cor.res
data=data[data$Name %in% SNP$annot,]
data$Zscore=data$Coefficient/data$Coefficient_std_error

data=data[data$Zscore>0,]
cor.res_temp=cor.res[rownames(cor.res) %in% data$Name, colnames(cor.res) %in% data$Name]
cor.res_temp=as.matrix(cor.res_temp)


data=data[match(rownames(cor.res_temp), data$Name), ]


if (nrow(data>1)) {

rss=susie_rss(R=cor.res_temp, z=data$Zscore, estimate_residual_variance = F, L=10, n=1190321)

pip=as.data.frame(rss$pip)
pip$Name=data$Name
pip$Zscore=data$Zscore
pip$Coeff=data$Coefficient
pip$TRAIT=trait
pip$INDEX=1:nrow(pip)
colnames(pip)[1]='PIP'
pip=pip[order(-pip[,1]),]
all_pips=rbind(all_pips,pip)





all_cs=data.frame()

a=susie_get_cs(rss,Xcorr=cor.res_temp,coverage=0.96,min_abs_corr = 0.3000922) # maybe use Xcorr to filter out based on correlation threshold
                    #Xcorr : p by p matrix of correlations between variables (covariates). When provided, it will be used to remove CSs whose minimum correlation among variables is smaller than min_abs_corr.

if (length(a$cs)==0) {print(paste0(rs,' no CS_95 with purity threshold!'));

} else for (j in 1:length(a$cs)) {
    cs=data.frame( 'CS_Num'=names(a$cs)[j],'Size'=length(a$cs[[j]]),'annots'=a$cs[[j]], 'TRAIT'=trait)
    cs$annot_name=pip$Name[match(cs$annots,pip$INDEX)]
    cs$PIP=pip$PIP[match(cs$annots,pip$INDEX)]
    cs$Zscore=pip$Zscore[match(cs$annots,pip$INDEX)]
    cs=cs[,-3]
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
    cs=cs[,c('CS_Num','Size','TRAIT','annot_name','PIP','Zscore','Type')]
    cs$S=S_all$S[match(cs$CS_Num,S_all$L)]
    cs_alpha=rbind(cs_alpha,cs)

    } else { print(paste(L$L[i],' already in Susie CS')) }  }
    } else {
        print('no alpha CS')
       }

print(paste0('Number of alpha CS: ', nrow(table(cs_alpha$CS_Num))))

all_cs=rbind(all_cs,cs_alpha)



## Write out

rs=basename(rs)
rs=gsub(".out","",rs)

outdir=paste0(wd,'/out/ctfmsnp_cancer/',trait)

print(paste0("TRAIT: ",trait," SNP: ",rs, " OUTDIR: ",outdir ))
print(paste0("number CS:",nrow(table(all_cs$CS_Num))))


system(paste0('mkdir -p ', outdir))

write.table(all_pips, file=paste0(outdir,'/',rs,'_pips.txt'), sep='\t', quote=F, row.names=F)
write.table(all_cs, file=paste0(outdir,'/',rs,'_CS.txt'), sep='\t', quote=F, row.names=F)


} else {print('Less than 1 annotation with postive Z')}