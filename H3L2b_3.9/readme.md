# generate qa code
- 1, generate minitree from Data and embedding (StRoot )
- 2, use runana.sh to read minitree and produce histograms
- 3, root -b drawData.C

# spectra calculation:
- 1, read minitree and generate histograms
- 2, root -b -q drawData_00_10.C
