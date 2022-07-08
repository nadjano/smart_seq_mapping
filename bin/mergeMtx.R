#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(DropletUtils))

cl <- commandArgs(trailingOnly = TRUE)

dirlist <- cl[1]
outdir <- cl[2]

if ( ! file.exists(dirlist)){
    stop(paste0('Directory listing file', dirlist, 'does not exist'))
}

dirs <- readLines(dirlist)

for ( dir in dirs){
    if ( ! dir.exists(dir)){
        stop(paste0('Provided directory', dir, 'does not exist'))
    }
}

mat <- read10xCounts(dirs[1])

if (length(dirs) > 1){

    for (i in 2:length(dirs)){
        print(paste0('Reading sub-directory ', i, '...'))

        dir <- dirs[i]
        nextmat <- read10xCounts(dir)
        intersectrows = intersect(rownames(mat), rownames(nextmat))

        mat <- cbind(mat[intersectrows,], nextmat[intersectrows,])
    }
}

print('Writing outputs...')
write10xCounts(outdir, assays(mat)[[1]], gene.id=rownames(mat), barcodes=colData(mat)$Barcode, overwrite=TRUE)
print('Done')