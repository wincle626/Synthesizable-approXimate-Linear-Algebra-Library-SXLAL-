# Synthesizable approXimate Linear Algebra Library ( SXLAL )
This is a C++ template libray for various precision linear algebra. It is aiming to research on the approximate computing for optimization algorithms with linear algebra arithmetics in the software layer. The code is synthesizable using High Lever Synthesis (HLS) tool targeting on the hardware design on FPGA. It's also extending to more Machine Learning applications. 

## Key Publications
[Efficient Reconfigurable Mixed Precision ℓ1 Solver for Compressive Depth Reconstruction. Journal of Signal Processing Systems 94, 1083–1099 (2022).](https://doi.org/10.1007/s11265-022-01766-3)

[Energy Efficient Approximate 3D Image Reconstruction. IEEE Transactions on Emerging Topics in Computing (PrePrints)](https://doi.ieeecomputersociety.org/10.1109/TETC.2021.3116471)

## Dependent Libraries
sudo apt install python python-numpy python-matplotlib libboost-all-dev libeigen3-dev fftw-dev -y

Notice that: currently the fixed-point arithmetics are all based on the Xilinx Vivado HLS fixed-point template. **You MUST have a Xilinx license first**. Then copy the "include" folder to the "header" folder and renamed as "xilinx". Otherwise, you cannot pass the compilation with fixed point support.

## Usage
### 1. Go to the build folder: cd build

### 2. Cmake the source folder: cmake ..

### 3. Build the source: make -j4

### 4. Excuete the example program: e.g. ./pgd

## Proximal Gradient Descend Example (PGD)
Different precision of gradient process is evaluated on x86 platform with Intel i7 6820HK with 8Gb DDR4 memory. (W: the total bit width, I: the integer bit width, M: the mantissa bit witdh)
### 1. Data size 64
<p align="center">
  float<img src="https://github.com/wincle626/SXLAL/blob/master/CMakeLib/data/x86_x64/PGD/64x64/float/ProximalGradientDecent.png" width="350">
  double<img src="https://github.com/wincle626/SXLAL/blob/master/CMakeLib/data/x86_x64/PGD/64x64/double/ProximalGradientDecent.png" width="350">
</p>

float: Number of iterations = 40529, Time Usage = 252.721 ms

double: Number of iterations = 40683, Time Usage = 574.483 ms

### 2. Data size 128
<p align="center">
  float<img src="https://github.com/wincle626/SXLAL/blob/master/CMakeLib/data/x86_x64/PGD/128x128/float/ProximalGradientDecent.png" "float" width="350">
  double<img src="https://github.com/wincle626/SXLAL/blob/master/CMakeLib/data/x86_x64/PGD/128x128/double/ProximalGradientDecent.png" "double" width="350">
</p>

float: Number of iterations = 47934, Time Usage = 833.31 ms

double: Number of iterations = 48191, Time Usage = 1487.99 ms

### 3. Data size 256
<p align="center">
  float<img src="https://github.com/wincle626/SXLAL/blob/master/CMakeLib/data/x86_x64/PGD/256x256/float/ProximalGradientDecent.png" width="350">
  double<img src="https://github.com/wincle626/SXLAL/blob/master/CMakeLib/data/x86_x64/PGD/256x256/double/ProximalGradientDecent.png" width="350">
</p>

float: Number of iterations = 224087, Time Usage = 11561.9 ms

double: Number of iterations = 224665, Time Usage = 21787.5 ms

### 4. Data size 512
<p align="center">
  float<img src="https://github.com/wincle626/SXLAL/blob/master/CMakeLib/data/x86_x64/PGD/512x512/float/ProximalGradientDecent.png" width="350">
  double<img src="https://github.com/wincle626/SXLAL/blob/master/CMakeLib/data/x86_x64/PGD/512x512/double/ProximalGradientDecent.png" width="350">
</p>

float: Number of iterations = 403136, Time Usage = 66061.3 ms

double: Number of iterations = 402571, Time Usage = 122316 ms

### 5. Data size 1024
<p align="center">
  float<img src="https://github.com/wincle626/SXLAL/blob/master/CMakeLib/data/x86_x64/PGD/1024x1024/float/ProximalGradientDecent.png" width="350">
  double<img src="https://github.com/wincle626/SXLAL/blob/master/CMakeLib/data/x86_x64/PGD/1024x1024/double/ProximalGradientDecent.png" width="350">
</p>

float: Number of iterations = 297133, Time Usage = 313763 ms

double: Number of iterations = 261360, Time Usage = 395742 ms

## Illustration of customized floating point precision agains data size

### 1. Data size 64

<p align="center">
  W10, M5<img src="https://github.com/wincle626/SXLAL/blob/master/CMakeLib/data/x86_x64/PGD/64x64/xfpt_w10_m5/ProximalGradientDecent_gfloat.png" width="350">
  W12, M5<img src="https://github.com/wincle626/SXLAL/blob/master/CMakeLib/data/x86_x64/PGD/64x64/xfpt_w12_m5/ProximalGradientDecent_gfloat.png" width="350">
  W16, M5<img src="https://github.com/wincle626/SXLAL/blob/master/CMakeLib/data/x86_x64/PGD/64x64/xfpt_w16_m5/ProximalGradientDecent_gfloat.png" width="350">
  W24, M9<img src="https://github.com/wincle626/SXLAL/blob/master/CMakeLib/data/x86_x64/PGD/64x64/xfpt_w24_m9/ProximalGradientDecent_gfloat.png" width="350">
</p>

W10, M5: Number of iterations = 99998, Time Usage = 11750 ms

W12, M5: Number of iterations = 99998, Time Usage = 11062.5 ms

W16, M5: Number of iterations = 99998, Time Usage = 11265.6 ms

W24, M9: Number of iterations = 40570, Time Usage = 4546.88 ms

### 2. Data size 128

<p align="center">
  W10, M5<img src="https://github.com/wincle626/SXLAL/blob/master/CMakeLib/data/x86_x64/PGD/128x128/xfpt_w10_m5/ProximalGradientDecent_gfloat.png" width="350">
  W12, M5<img src="https://github.com/wincle626/SXLAL/blob/master/CMakeLib/data/x86_x64/PGD/128x128/xfpt_w12_m5/ProximalGradientDecent_gfloat.png" width="350">
  W16, M5<img src="https://github.com/wincle626/SXLAL/blob/master/CMakeLib/data/x86_x64/PGD/128x128/xfpt_w16_m5/ProximalGradientDecent_gfloat.png" width="350">
  W24, M9<img src="https://github.com/wincle626/SXLAL/blob/master/CMakeLib/data/x86_x64/PGD/128x128/xfpt_w24_m9/ProximalGradientDecent_gfloat.png" width="350">
</p>

W10, M5: Number of iterations = 99998, Time Usage = 41703.1 ms

W12, M5: Number of iterations = 99998, Time Usage = 41781.2 ms

W16, M5: Number of iterations = 99998, Time Usage = 35000 ms

W24, M9: Number of iterations = 47960, Time Usage = 19796.9 ms

## Illustration of fixed-point precision agains data size

### 1. Data size 64

<p align="center">
  W28, I10<img src="https://github.com/wincle626/SXLAL/blob/master/CMakeLib/data/x86_x64/PGD/64x64/fxpt_w28_i10/ProximalGradientDecent.png" width="350">
  W32, I12<img src="https://github.com/wincle626/SXLAL/blob/master/CMakeLib/data/x86_x64/PGD/64x64/fxpt_w32_i12/ProximalGradientDecent.png" width="350">
  W48, I12<img src="https://github.com/wincle626/SXLAL/blob/master/CMakeLib/data/x86_x64/PGD/64x64/fxpt_w48_i12/ProximalGradientDecent.png" width="350">
</p>

W28, I10: Number of iterations = 1617, Time Usage = 597.451 ms

W32, I12: Number of iterations = 3678, Time Usage = 3830.12 ms

W48, I12: Number of iterations = 41602, Time Usage = 70051.8 ms

### 2. Data size 128

<p align="center">
  W28, I12<img src="https://github.com/wincle626/SXLAL/blob/master/CMakeLib/data/x86_x64/PGD/128x128/fxpt_w28_i12/ProximalGradientDecent.png" width="350">
  W32, I12<img src="https://github.com/wincle626/SXLAL/blob/master/CMakeLib/data/x86_x64/PGD/128x128/fxpt_w32_i12/ProximalGradientDecent.png" width="350">
  W48, I16<img src="https://github.com/wincle626/SXLAL/blob/master/CMakeLib/data/x86_x64/PGD/128x128/fxpt_w48_i16/ProximalGradientDecent.png" width="350">
</p>

W28, I12: Number of iterations = 2635, Time Usage = 332.926 ms

W32, I12: Number of iterations = 3228, Time Usage = 877.335 ms

W48, I16: Number of iterations = 36661, Time Usage = 17046.1 ms

## Various precision against data size

As shown, the iteration numbers using float and double precision across different dimention size are very similar from 64 to 256. When reaching the dimension 1024, there is a siginificant difference of iteration number between float and double precision data type. 
<p align="center">
  <img src="https://github.com/wincle626/SXLAL/blob/master/CMakeLib/data/figs/Converge_Iter_Compare_fvd.png" width="350">
</p>

As shown, the executing time of double precision is alway larger than the float precision even when less iteration happens at dimension 1024. 
<p align="center">
  <img src="https://github.com/wincle626/SXLAL/blob/master/CMakeLib/data/figs/Converge_Time_Compare_fvd.png" width="350">
</p>

## Various precision under the same data size

### 1. Data size 64

As shown, the iteration number of fixed point precision is significantly lower than the floating point precision but it is at a cost of degree of optimization. To be specifically, the gradient descend iteration termiated earlier due to the quantization of the either integer or the fraction part of the fixed point number. It is hard to garantee there is both enough bit to represent the integer or the fraction part unless very large number of bits are used. So when quantized, the gradient decent cannot continue any more as the error between current iteration and last iteration will be ZERO, which means there is no further 'descend'. 

The executing time comparison between fixed point and floating point is not fair enough on x86 hard core processor as the fixed point has large overhead in approximating the floating number into fixed point. The fixed point arithmetics are based on the Xilinx HLS library headers. The results are better to be considered as a verification of the various precision effect on the proximal gradient descend algorithm. 

By comparing the optimized solution between C++ and CVX code, it is clear that both float and double precision maitains the better optimized solution to the problem with only 0.0007 normalized error compared to CVX solution, while the fixed point precision has much wrose normalized error up to 0.25 compared to CVX solution. By using the fixed point precision, minimun value of cost function is much larger than the float and double precision due to the fact that it has not reash as optimal as float and double precision did when doing the gradient descend. 

<p align="center">
  <img src="https://github.com/wincle626/SXLAL/blob/master/CMakeLib/data/figs/Converge_Iter_Compare_64x64.png" width="350">
  <img src="https://github.com/wincle626/SXLAL/blob/master/CMakeLib/data/figs/Converge_Time_Compare_64x64.png" width="350">
  <img src="https://github.com/wincle626/SXLAL/blob/master/CMakeLib/data/figs/Opt_Err_Compare_64x64.png" width="350">
  <img src="https://github.com/wincle626/SXLAL/blob/master/CMakeLib/data/figs/Opt_Solu_Compare_64x64.png" width="350">
</p>

### 2. Data size 128

Similar situation happens for the dimension 128. 
<p align="center">
  <img src="https://github.com/wincle626/SXLAL/blob/master/CMakeLib/data/figs/Converge_Iter_Compare_128x128.png" width="350">
  <img src="https://github.com/wincle626/SXLAL/blob/master/CMakeLib/data/figs/Converge_Time_Compare_128x128.png" width="350">
  <img src="https://github.com/wincle626/SXLAL/blob/master/CMakeLib/data/figs/Opt_Err_Compare_128x128.png" width="350">
  <img src="https://github.com/wincle626/SXLAL/blob/master/CMakeLib/data/figs/Opt_Solu_Compare_128x128.png" width="350">
</p>

## Alternative Directive Method of Multiplier (ADMM)

### Eclipse Time of 4000 iterations

#### 1.Eclipse Time of Matlab: 8.937781 seconds

#### 2. Eclipse Time of C/C++: 101.632 seconds (float), 102.163 seconds (double)

### Relative Error v.s. CVX

#### 1. Matlab: 0.040783

#### 2. C/C++: 0.040775 (float), 0.040783 (double)

#### This matches the benchmark of speed comparison between C/C++ and BLAS implementation in fortran at [Eigen3 official website](https://eigen.tuxfamily.org/index.php?title=Benchmark). 

#### Matlab is using BLAS since Day 1, which is why the linear algebra, like matrix operations, is much faster than C/C++ implementation. 

## Some Inspects

### 1. The double and float precision provide fairly good accurate solution compared to the CVX solution.

### 2. As the dimension grows, it is hard to use the fixed point precision as the sum of values increases the requirement of integer part bitwidth while relatively large fraction part bitwidth is still required to maintain the gradient descend steps. 

### 3. The fixed point precision is not suitable to run barely on hard core processor as most of the modern processor are floating point precision and does not support fixed point instructions. This will cause a lot overhead to emulate the fixed point arithmetic operators. 

### 4. Floating point still has its advantage in representing the number in a relative large dynmaic value range. But there is no arbitray floating point format exist on modern processor. 

### 5. There is universal numbers, UNUM, data format alternative to the IEEE 754 standard. It worthes a exploration but again most of the modern processor does not support UNUM arithmetic instructions. It is doubt that the same situation happens to it as the fixed point data type. We will see. 

### 6. Through the ADMM example, it is very clear that in order to gain speedup of linear algebra, the simple approximation, such as various precision is not enough. Approximation does reduce the power consumption but it might also degrade the performance of either algorithm or executing time. So a concurrrent/parallel/distributed architecture is a must in such cases that reduce the executing time and improve the algorithm performance by adopting more executions, even with approximated approach. 

### 7. The customized floating point precision is not supported by any modern processor which means the meaning of it is the same as the fixed point precision. It can ONLY be used as simulator for evaluating the precision effect on the algorithm but not consider the energy effect in practice. 

## Reference

### 1. [Approximate Computing: Challenges And Opportunities. ](https://ieeexplore.ieee.org/document/7738674)

### 2. [Approximate Computing: A Survey.](https://ieeexplore.ieee.org/document/7348659)

### 3. [A Survey Of Techniques for Approximate Computing.](https://dl.acm.org/citation.cfm?id=2893356)

## External Links

### 1. [Xilinx HLS user guild.](https://www.xilinx.com/support/documentation/sw_manuals/xilinx2017_4/ug902-vivado-high-level-synthesis.pdf)

### 2. [Xilinx fixed point white paper.](https://www.xilinx.com/support/documentation/white_papers/wp491-floating-to-fixed-point.pdf)

### 3. [CVX: Matlab Software for Disciplined Convex Programming.](http://cvxr.com/cvx/)

### 4. [OPRECOMP project: aiming to build an innovative, reliable foundation for computing based on transprecision analytics.](http://oprecomp.eu/)

### 5. [Eigen3: a C++ template library for linear algebra: matrices, vectors, numerical solvers, and related algorithms.](http://eigen.tuxfamily.org/index.php?title=Main_Page)

### 6. [Posit & Unum: next generation arithmetic](https://posithub.org/index)
