# Description

An exhaustive method based on high-performance GPU computing to reveal degenerate oligonucleotide motifs written in the 15-letter IUPAC code.

# Build 
```
git clone git@github.com:friofry/motif_finder_cuda.git
cd motif_finder_cuda
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
cmake --build . -j
benchmark_cpu/benchmark_cpu 
```

# Run tests
```
cd ../test_data/test_0
../../build/motif_finder_cpu/motif_finder_cpu
```