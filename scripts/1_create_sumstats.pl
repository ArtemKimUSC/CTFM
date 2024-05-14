use Math::CDF; #conda install -c bioconda perl-math-cdf
use Getopt::Long;
use FindBin qw($RealBin);

$beta_col=""; $se_col=""; $Z_col=""; $pval_col=""; $OR_col=""; $myN=0; $N_col="";

GetOptions ("filein=s" => \$FILEIN, "fileout=s" => \$FILEOUT, "chr=i" => \$chr_col, "pos=i" => \$pos_col, "rs=i" => \$rs_col, "A1=i" => \$A1_col, "A2=i" => \$A2_col, "beta=i" => \$beta_col, "se=i" => \$se_col, "Z=i" => \$Z_col, "pval=i" => \$pval_col, "OR=i" => \$OR_col, "myN=i" => \$myN, "N=i" => \$N_col, "N1=i" => \$N1_col, "N2=i" => \$N2_col, "hg=i" => \$hg);

printf("perl create_sumstats.pl\n");
printf("     --filein $FILEIN\n");
printf("     --fileout $FILEOUT\n");
if ($rs_col==""){
printf("     --chr $chr_col\n");
printf("     --pos $pos_col\n");
} else {
printf("     --rs $rs_col\n");
}
printf("     --A1 $A1_col\n");
printf("     --A2 $A2_col\n");
#
if ($beta_col!=""){
printf("     --beta $beta_col\n");
printf("     --se $se_col\n");
} elsif ($Z_col!=""){
printf("     --Z $Z_col\n");
} else {
printf("     --pval $pval_col\n");
if ($OR_col!=""){
printf("     --OR $OR_col\n");
}
}
#
if ($myN!=0){
printf("     --myN $myN\n");
} elsif ($N_col!=""){
printf("     --N $N_col\n");
} else {
printf("     --N1 $N1_col\n");
printf("     --N2 $N2_col\n");
}
printf("     --hg $hg\n");

##1) Check SNPs in bim + weight
print "\nRead HM3 SNPs\n";
#weight
open(IN,"$RealBin/../data/w_hm3.txt" ) || die ;
while (<IN>) {
    chomp $_; @line=split;
    $rsin{$line[0]}=1;
}
close IN;
#bim
if ($hg==19){$ref_file="$RealBin/../data/1000G_EUR_Phase3/plink_files/1000G.EUR.QC"}
if ($hg==38){$ref_file="$RealBin/../data/1000G_EUR_Phase3_hg38/plink_files/1000G.EUR.hg38"}
print "Read SNP info in bim files\n";
for (my $chr = 1; $chr <= 22; $chr++) { 
    open(IN,"$ref_file.$chr.bim" ) || die ;
	while (<IN>) {
		chomp $_; @line=split;
        $alleles = "$line[4]$line[5]";
        if (($alleles ne "AT") && ($alleles ne "TA") && ($alleles ne "CG") && ($alleles ne "GC")) {
            if (defined($rsin{$line[1]})){
                if ($rs_col==""){
                    $rs{$chr}{$line[3]}{$line[4]}{$line[5]}=$line[1];
                    $rs{$chr}{$line[3]}{$line[5]}{$line[4]}=$line[1];
                } else {
                    $rs{$line[1]} = 1;
                }
            }
        }
	}
	close IN;
}

#2) Create new my sumstats file
print "\nCreate $FILEOUT \n";
open(OUT,">$FILEOUT");
printf OUT "SNP	A1	A2	N	Z\n";
system("gunzip $FILEIN.gz");
open(IN,"$FILEIN" ) || die ;
	while (<IN>) {
    chomp $_; @line=split;
    if ($rs_col==""){
        if ( defined($rs{$line[$chr_col-1]}{$line[$pos_col-1]}{uc($line[$A1_col-1])}{uc($line[$A2_col-1])})){
            $thisSNP = $rs{$line[$chr_col-1]}{$line[$pos_col-1]}{uc($line[$A1_col-1])}{uc($line[$A2_col-1])};
            $thisA1 =  uc($line[$A1_col-1]);
            $thisA2 =  uc($line[$A2_col-1]);
            if ($myN!=0){
                $thisN = $myN
            } elsif ($N_col!=""){
                $thisN = $line[$N_col-1]
            } else {
                $thisN = $line[$N1_col-1]+$line[$N2_col-1]
            }
            if ($beta_col!=""){
                if ($line[$se_col-1]!=0){
                    $thisZ = sprintf("%.4f", $line[$beta_col-1]/$line[$se_col-1]);
                }
            } elsif ($Z_col!=""){
                $thisZ = sprintf("%.4f", $line[$Z_col-1]);
            } else {
                $thisZ = sprintf("%.4f", &Math::CDF::qnorm($line[$pval_col-1]/2)*-1);
                if ($OR_col!=""){ if ($line[$OR_col-1]<1){  $thisZ = -1 *$thisZ } }
            }
            if (($beta_col!="") && ($line[$se_col-1]!=0)){
                printf OUT "$thisSNP	$thisA1	$thisA2	$thisN	$thisZ\n";    
            }
        }
    } else {
        if ( defined($rs{$line[$rs_col-1]})){
            $thisSNP = $line[$rs_col-1];
            $thisA1 =  uc($line[$A1_col-1]);
            $thisA2 =  uc($line[$A2_col-1]);
            if ($myN!=0){
                $thisN = $myN
            } elsif ($N_col!=""){
                $thisN = $line[$N_col-1]
            } else {
                $thisN = $line[$N1_col-1]+$line[$N2_col-1]
            }
            if ($beta_col!=""){
                if ($line[$se_col-1]!=0){
                    $thisZ = sprintf("%.4f", $line[$beta_col-1]/$line[$se_col-1]);
                }
            } elsif ($Z_col!=""){
                $thisZ = sprintf("%.4f", $line[$Z_col-1]);
            } else {
                $thisZ = sprintf("%.4f", &Math::CDF::qnorm($line[$pval_col-1]/2)*-1);
                if ($OR_col!=""){ if ($line[$OR_col-1]<1){  $thisZ = -1 *$thisZ } }
            }
            #if (($beta_col!="") && ($line[$se_col-1]!=0)){
                printf OUT "$thisSNP	$thisA1	$thisA2	$thisN	$thisZ\n";    
            #}
        }
    } 
}
close IN;
close OUT;
system("gzip $FILEIN");
system("gzip $FILEOUT");
