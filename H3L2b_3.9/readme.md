# generate qa code
- 1, generate minitree from Data and embedding (StRoot), RCF address: /star/u/jiyj/pwg/Hypertron/data/7p3_2020_2body
- 2, use runana.sh to read minitree and produce histograms
- 3, root -b drawData.C

# spectra calculation:
- 1, read minitree and generate histograms
- 2, root -b -q drawData_00_10.C
- 3, Update on May 06, add TGraphErrors* ptshift(TF1* f, TH1* h) function, f is the function to fit the spectra and h is the original spectra data with orignal pT binning.
    Th calculation is the same as https://www.star.bnl.gov/protected/lfspectra/marr/Analysis/PlaceDataPoints.pdf
