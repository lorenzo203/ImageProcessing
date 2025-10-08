# NLA Challenge 1

This repository contains the source code for **Challenge 1**, which uses the [LIS library](https://www.ssisc.org/lis/) and [Eigen3](https://eigen.tuxfamily.org/) to solve sparse linear systems and produce output images.

## Requirements

Before building the code, make sure the following are available on your system:

- A C++17 compiler (e.g. `g++`)
- `make`
- [Eigen3](https://eigen.tuxfamily.org/) headers installed system-wide  
  (for example, on Debian/Ubuntu: `sudo apt install libeigen3-dev`)

The LIS library source (`lis-2.1.10/`) is already included in this repository.  
There is no need to install it manually: the `Makefile` takes care of all required LIS build commands automatically.

## Building

From the project root directory, run:

```bash
make
```

This will:

1. Configure and compile the LIS library (inside `lis-2.1.10/`) with the required options.

2. Compile `challenge1.cpp` with the appropriate include and library paths.

3. Link against LIS and Eigen.

4. Produce the executable `./challenge1`.

So `make` handles both the LIS build process and the compilation of the program in a single step.

## Running

After building, simply run:

```bash
./challenge1
```

## Cleaning up

To remove all compiled files and reset the directory to a clean state:

```bash
make clean
```

This setup ensures that anyone with Eigen installed can compile and run the project using just `make` and `./challenge1`, without needing to install LIS separately.

## Repository Structure

```
thisrepo/  
├── challenge1.cpp # Main C++ source code  
├── stb_image_write.h # Header for writing images (PNG)  
├── stb_image.h # Header for reading images  
├── Outputs/ # Folder where generated outputs are saved  
├── lis-2/ # LIS library source files used for compilation  
├── Makefile # Build script that compiles LIS and challenge1.cpp  
└── uma.jpg # Input image used by the program
```

## Output

The output obtained by the execution of this program is:

```
Task 1. Image size: 500 x 456 (rows x cols)
Task 2. Noise image saved to Outputs/noise.png
Task 3. v size = 228000, ||v||_2 = 62506.566799
Task 4. Non-zero entries in A1 = 1592178
Task 5. Smoothed image saved to Outputs/smoothing.png
Task 6. Non-zero entries in A2 = 1138088, A2 symmetric? true
Task 7. Sharpened image saved to Outputs/sharpening.png
Task 8. Exported A2.mtx and w.mtx under /Outputs
Task 8. LIS iterations = 11, final residual = 9.35333117183e-16
Task 9. LIS solution image saved to Outputs/lis_solution.png
Task 10. A3 symmetric? false
Task 11. Edge detection image saved to Outputs/edge_detection.png
Task 12. Eigen solver iterations = 57, final residual = 4.98074833471e-09
Task 13. Eigen solution image saved to Outputs/eigen_solution.png
```
