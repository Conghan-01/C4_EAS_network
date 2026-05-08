
DEBUG <- FALSE;
VERBOSE <- FALSE;
MAX_QUAL <- 99.0;

run <- function(inputFile, category="ALL") {
    if (VERBOSE) {
        cat(sprintf("%s: Reading input data ...\n", date()));
    }
    inputData = read.table(inputFile, header=T, sep="\t", stringsAsFactors=F, check.names=F, row.names=1);
    pcnData = apply(inputData[2:nrow(inputData),,drop=F], 1:2, function(v) { strsplit(v,";",fixed=T)[[1]][[1]] });
    pcnlData = apply(inputData[2:nrow(inputData),,drop=F], 1:2, function(v) { strsplit(v,";",fixed=T)[[1]][[2]] });
    if (category != "ALL") {
        pcnData = pcnData[,inputData["PST",] == category];
        pcnlData = pcnlData[,inputData["PST",] == category];
        if (ncol(pcnData) == 0) {
            cat(sprintf("No data for PST category %s\n", category));
            return(invisible());
        }
    }

    if (VERBOSE) {
        cat(sprintf("%s: Refining haplotypes ...\n", date()));
    }
    results = NULL;
    samples = colnames(pcnData);
    for (sample in samples) {
        for (hap in 1:2) {
            hapData = processHaplotype(sample, hap, pcnData, pcnlData);
            results = rbind(results, hapData);
        }
    }

    outData = results;
    outData$OQUAL = sprintf("%1.2f", outData$OQUAL);
    outData$RQUAL = sprintf("%1.2f", outData$RQUAL);
    write.table(outData, quote=F, sep="\t", row.names=F);
    if (VERBOSE) {
        cat(sprintf("%s: Analysis done.\n", date()));
    }
}

processHaplotype <- function(sample, hap, pcnData, pcnlData) {
    if (DEBUG) {
        cat(sprintf("#DBG: %s processing %s-%d ...\n", date(), sample, hap));
    }
    cn_C4 = getCN(pcnData, "CNV_C4", sample, hap);
    cn_C4A = getCN(pcnData, "CNV_C4A", sample, hap);
    cn_C4B = getCN(pcnData, "CNV_C4B", sample, hap);
    cn_HERV = getCN(pcnData, "CNV_HERV", sample, hap);
    cnls_C4 = getCNLs(pcnlData, "CNV_C4", sample, hap);
    cnls_C4A = getCNLs(pcnlData, "CNV_C4A", sample, hap);
    cnls_C4B = getCNLs(pcnlData, "CNV_C4B", sample, hap);
    cnls_HERV = getCNLs(pcnlData, "CNV_HERV", sample, hap);
    originalCNs = c(cn_C4, cn_C4A, cn_C4B, cn_HERV);
    refinedHaplotype = refineHaplotype(cnls_C4, cnls_C4A, cnls_C4B, cnls_HERV);
    refinedCNs = refinedHaplotype$CNS;
    originalHap = sprintf("H_%d_%d_%d_%d", originalCNs[[1]], originalCNs[[2]], originalCNs[[3]], originalCNs[[4]]);
    refinedHap = sprintf("H_%d_%d_%d_%d", refinedCNs[[1]], refinedCNs[[2]], refinedCNs[[3]], refinedCNs[[4]]);
    originalQual = computeQual(originalCNs, cnls_C4, cnls_C4A, cnls_C4B, cnls_HERV);
    refinedQual = refinedHaplotype$QUAL;
    result = data.frame(SAMPLE=sample, HAP=hap, ORIGINAL=originalHap, REFINED=refinedHap, OQUAL=originalQual, RQUAL=refinedQual);
    return(result);
}

refineHaplotype <- function(cnls_C4, cnls_C4A, cnls_C4B, cnls_HERV) {
    hmax1 = NULL;
    llmax1 = -Inf;
    llmax2 = -Inf;
    N = length(cnls_C4);
    for (cn_C4 in seq(0, N-1)) {
        for (cn_C4A in seq(0, cn_C4)) {
            for (cn_HERV in seq(0, cn_C4)) {
                cn_C4B = cn_C4 - cn_C4A;
                ll_C4 = getCNL(cnls_C4, cn_C4);
                ll_C4A = getCNL(cnls_C4A, cn_C4A);
                ll_C4B = getCNL(cnls_C4B, cn_C4B);
                ll_HERV = getCNL(cnls_HERV, cn_HERV);
                ll_hap = ll_C4 + ll_C4A + ll_C4B + ll_HERV;
                if (ll_hap >= llmax1) {
                    llmax2 = llmax1;
                    llmax1 = ll_hap;
                    hmax1 = c(cn_C4, cn_C4A, cn_C4B, cn_HERV);
                } else if (ll_hap >= llmax2) {
                    llmax2 = ll_hap;
                }
            }
        }
    }
    if (is.infinite(llmax1) || is.infinite(llmax2)) {
        qual = 0;
    } else {
        qual = pmin(10 * (llmax1 - llmax2), MAX_QUAL);
    }
    result = list(CNS=hmax1, QUAL=qual);
    return(result);
}

computeQual <- function(cns, cnls_C4, cnls_C4A, cnls_C4B, cnls_HERV) {
    llmax1 = -Inf;
    llmax2 = -Inf;
    N1 = length(cnls_C4);
    N2 = length(cnls_C4A);
    N3 = length(cnls_C4B);
    N4 = length(cnls_HERV);
    for (cn_C4 in seq(0, N1-1)) {
        for (cn_C4A in seq(0, N2-1)) {
            for (cn_C4B in seq(0, N3-1)) {
                for (cn_HERV in seq(0, N4-1)) {
                    ll_C4 = getCNL(cnls_C4, cn_C4);
                    ll_C4A = getCNL(cnls_C4A, cn_C4A);
                    ll_C4B = getCNL(cnls_C4B, cn_C4B);
                    ll_HERV = getCNL(cnls_HERV, cn_HERV);
                    ll = ll_C4 + ll_C4A + ll_C4B + ll_HERV;
                    if (all(c(cn_C4, cn_C4A, cn_C4B, cn_HERV) == cns)) {
                        llmax1 = ll;
                    } else if (ll > llmax2) {
                        llmax2 = ll;
                    }
                }
            }
        }
    }
    if (is.infinite(llmax1) || is.infinite(llmax2) || llmax1 < llmax2) {
        qual = 0;
    } else {
        qual = pmin(10 * (llmax1 - llmax2), MAX_QUAL);
    }
    return(qual);
}

getCN <- function(pcnTable, site, sample, hap) {
    pcnText = pcnTable[site,sample];
    pcnText = strsplit(pcnText, "|", fixed=T)[[1]][[hap]]
    return(as.numeric(pcnText))
}

getCNLs <- function(pcnlTable, site, sample, hap) {
    pcnlText = pcnlTable[site,sample];
    pcnlText = strsplit(pcnlText, "|", fixed=T)[[1]][[hap]]
    return(as.numeric(strsplit(pcnlText, ",", fixed=T)[[1]]))
}

getCNL <- function(cnlvec, cn) {
    if (cn+1 > length(cnlvec)) {
        return(0);
    }
    return(cnlvec[cn+1]);
}

args = commandArgs(TRUE)
do.call(run, as.list(args))
