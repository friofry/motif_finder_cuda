0. Download cmake 3.20 (already done):
   * https://github.com/Kitware/CMake/releases/download/v3.20.3/cmake-3.20.3-linux-x86_64.tar.gz
   * ```tar -xvzf cmake-3.20.3-linux-x86_64.tar.gz```
    
1. Download and unpack source files:
https://github.com/friofry/motif_finder_cuda/archive/refs/tags/v1.0.zip


2. Run the a6500g10q node:
```shell
qsub -I -l walltime=00:30:00 -l select=1:ngpus=4:ncpus=12:mem=40gb -q a6500g10q@vm-pbs2
```
3. Unpack the sources:
```shell
unzip v1.0.zip -d .
```
4. Navigate to the sourcedir:
```shell
cd motif_finder_cuda
```
5. Create build dir and navigate to it:
```shell
mkdir biuld
cd build
```
6. Configure project:
```shell
~/cmake-3.20.3-linux-x86_64/bin/cmake -DCMAKE_BUILD_TYPE=Release ..
```

7. Build:
```shell
/cmake-3.20.3-linux-x86_64/bin/cmake --build . -j
```

8. Run:
* Copy test data:
```asm
cp ../motif_finder_gpu/init.ini ../motif_finder_gpu/test_12.fst motif_finder_gpu/
```
* Run:
```asm
cd motif_finder_gpu
./motif_finder_gpu
```
