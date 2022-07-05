#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(DropletUtils))

commands = commandArgs(trailingOnly=T)

input = commands[1]
output_path = commands[2]

data = readRDS(input)
write10xCounts(x=data@assays$RNA@data, path=output_path)
write.table(data@meta.data, file.path(output_path, 'metadata'), sep='\t', quote=F)
