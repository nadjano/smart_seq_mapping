#!/usr/bin/env python

import argparse
import pandas as pd


parser = argparse.ArgumentParser(description= 'Produce cell medadata file from barcodes.tsv sdrf.file and cells.file.')
parser.add_argument('barcodes_file', help = 'barcodes.tsv file')
parser.add_argument('sdrfFile', help = 'sdrfFile')
parser.add_argument('cells_txt', help = 'file with cell_type information per cell')
parser.add_argument('out_file', help = 'name of output')
args = parser.parse_args() 

barcodes_file =args.barcodes_file
sdrfFile=args.sdrfFile
cells_txt=args.cells_txt
out = args.out_file

#read in files
barcodes = pd.read_csv(barcodes_file, header=None, names = ["Cell ID"])
sdrf = pd.read_csv(sdrfFile,header=0,  delimiter="\t")
cells = pd.read_csv(cells_txt, header=0,  delimiter="\t")

print(sdrf)
#seperate run ID and barcode
f = lambda x: x[0].split("-")[0]
barcodes["Comment[BioSD_SAMPLE]"] = barcodes.apply(f, axis=1)


#select relevant fields from sdrf  and cells file
sdrf_short = pd.DataFrame(sdrf[["Comment[BioSD_SAMPLE]", "Characteristics[individual]"]]).drop_duplicates()
cells = pd.DataFrame(cells[["Cell ID", "Inferred cell type - ontology labels"]])

#merge matching runs to get batch
merge = pd.merge(barcodes, sdrf_short, how = 'left',on="Comment[BioSD_SAMPLE]").drop_duplicates()

#merge matching Cell ID to get cell types
cell_meta = pd.merge(merge, cells, how = 'left',on="Cell ID").drop_duplicates(subset = "Cell ID")

#rename columns
cell_meta.columns = ["id", "run", "individual", "inferred_cell_type_-_ontology_labels"]

#write file
cell_meta.to_csv(out, sep='\t', index = False)