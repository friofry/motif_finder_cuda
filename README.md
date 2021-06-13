# Description

An exhaustive method based on high-performance GPU computing to reveal degenerate oligonucleotide motifs written in the 15-letter IUPAC code.

# Building from source code
## Dependencies
* CUDA 10 and higher (https://docs.nvidia.com/cuda/cuda-installation-guide-linux/index.html)
* Cmake 3.17  and higher (https://cmake.org/download/)
* GCC 6 and higher 

## Download
```
git clone git@github.com:friofry/motif_finder_cuda.git
```
or download the latest zip archive:
```
curl -L https://github.com/friofry/motif_finder_cuda/archive/refs/tags/v1.0.zip --output v1.0.zip
```

## Build 
Run these commands in the source code directory:
1. Make a build directory
```
mkdir build
cd build
```
2. Configure the project
```
cmake -DCMAKE_BUILD_TYPE=Release ..
```
To provide the exact path to NVIDIA CUDA SDK use the following command:
```
cmake -DCUDA_SDK_ROOT_DIR=<cuda_sdk_path> -DCMAKE_BUILD_TYPE=Release ..
```

3. Build 
```
cmake --build . -j
```

## motif_finder_cuda binary 

The binary is in `motif_finder_gpu` folder


**init.ini** file example:
```
1		Complementarity. 0 - forward strand. 1 - forward + reverse strand.
8		Minimum score
35		Maximum presence of motif for random reasons in the positive set of sequences [0-100]
5		Minimum presence of motif in the positive set of sequences [0-100]
test_12.fst	File with positive set of sequences
0		0 - neutral frequencies, 1 - real nucleotide frequencies in the set of sequences [0, 1]
0		[deprecated]
10		Maximum number of result motifs. 0 - reveal all significant motifs [0, ]
0		Markov chain order (0-Bernulli, 1-dinucleotide, 2-trinucleotide), when using real nucleotide frequencies [0-3]
0		[deprecated]
0		[deprecated]
0.fst		File with contrast set of sequences
90		Maximum score in a contrast set of sequences
0               [deprecated]
0               Output results with Bonferroni correction. 0 - without correction , 1 - use correction [0, 1]
a.txt           Output file
0               Output results in integer format (0 - real values, 1 - integer values). [0, 1]
1               Skip motifs that coincide with the previous ones with a shift. [0, 1]
```

## Run motif_finder_gpu

1. Put init.init, motif_finder_gpu and positive, negative sequences to the same directory
2. run `./motif_finder_gpu`

## Test data
1. Navigate to test folder:
```
cd ../test_data/test_0
```
2. Copy motif_finder_gpu to that folder
3. Run the motif_finder_gpu
```
./motif_finder_gpu
```

# Benchmarking the algorithm
See details [here](bencharking_performance.md)

# Running on NUSC cluster
See details [here](running_on_nusc_cluster.md)