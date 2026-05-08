# SCZ Heritability Analysis using MAGMA

zcat /gpfs/hpc/home/chenchao/hanc/project/03Cross_devlopment_eQTL/cross_develop/05intergration/GWAS/SCZ/pgc_scz_eas_autosome.maf.gz | awk '{
    if(NR==1) {
        print "SNP", "P", "N"
    } else {
        print $2, $11, $17+$18
    }
}' > pgc_scz_eas_magma.pval


zcat /gpfs/hpc/home/chenchao/hanc/project/03Cross_devlopment_eQTL/cross_develop/05intergration/GWAS/BD_2025raw/bip2024_eas_no23andMe.gz | awk '{
    if(NR==1) {
        print "SNP", "P", "N"
    } else {
        print $1, $9, $15+$16
    }
}' > bip_eas_magma.pval


zcat /gpfs/hpc/home/chenchao/hanc/project/03Cross_devlopment_eQTL/cross_develop/05intergration/GWAS/SCZ/pgc_scz_eas_autosome.maf.gz | awk '{
    if(NR==1) {
        print "SNP", "CHR", "BP"
    } else {
        print $2, $1, $3
    }
}' > pgc_scz_eas_magma.snp.loc


zcat /gpfs/hpc/home/chenchao/hanc/project/03Cross_devlopment_eQTL/cross_develop/05intergration/GWAS/BD_2025raw/bip2024_eas_no23andMe.gz| awk '{
    if(NR==1) {
        print "SNP", "CHR", "BP"
    } else {
        print $1, $2, $3
    }
}' > bip_eas_magma.snp.loc
