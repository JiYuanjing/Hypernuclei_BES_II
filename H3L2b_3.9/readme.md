# generate qa code
- 1, generate minitree from Data and embedding (StRoot), RCF address: /star/u/jiyj/pwg/Hypertron/data/7p3_2020_2body
- 2, use runana.sh to read minitree and produce histograms
- 3, root -b drawData.C

# spectra calculation:
- 1, read minitree and generate histograms
- 2, root -b -q drawData_00_10.C
