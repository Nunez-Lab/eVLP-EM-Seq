#!/usr/bin/env Rscript

library(DSS)
require(bsseq)

input_path = file.path("data", "methylation-dss-combined-filtered")

# ORIGINAL TOP STRAND

print("Reading S1...")
s1 = read.table(file.path(input_path, "WT_PJC_S1.txt"), header=TRUE)

print("Reading S2...")
s2 = read.table(file.path(input_path, "CD55off_PJC_S2.txt"), header=TRUE)

print("Reading S3...")
s3 = read.table(file.path(input_path, "WT_JPL_S3.txt"), header=TRUE)

print("Reading S4...")
s4 = read.table(file.path(input_path, "CD55off_JPL_S4.txt"), header=TRUE)

print("Making BSobj...")

BSobj = makeBSseqData(
    list(s2, s4, s1, s3),
    c("C2", "C4", "W1", "W3")
)

print("Running differential methylation test...")
dmlTest = DMLtest(
    BSobj,
    group1=c("C2", "C4"),
    group2=c("W1", "W3"),
    smoothing=TRUE,
    ncores=6
)

print("Writing CSV...")
write.csv(dmlTest, "data/dss-output/dss.csv", row.names=FALSE)
