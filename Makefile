NAME_FLAG =-o lbgpu
CUDA_ARCH_FLAG =-arch=sm_13

GCC_FLAGS =
GCC_FLAGS+=$(NAME_FLAG)
GCC_FLAGS+=-x c
GCC_FLAGS+=-g
GCC_FLAGS+=-p

all:
	gcc $(GCC_FLAGS) -D__device__='' -D__constant__='' ./src/lbgpu_prime.cu -lm

cuda:
	nvcc -v $(NAME_FLAG) ./src/lbgpu_prime.cu -lm $(CUDA_ARCH_FLAG)

mpi:
	mpicc -DPARALLEL=1 -o lb3d ./src/lb3d_prime.c -lm

warn:
	gcc -pedantic -Wall -Wunused-variable -o lb3d ./src/lb3d_prime.c

nocygwin:
	gcc -mno-cygwin -o lb3d ./src/lb3d_prime.c

sweep:
	/bin/rm -f ./out/*.dat
	/bin/rm -f ./out/rho*.txt
	/bin/rm -f ./out/u*.txt
	/bin/rm -f ./out/force*.txt
	/bin/rm -f ./out/f*.txt
	/bin/rm -f ./out/*.txt
	/bin/rm -f ./out/rho*.bmp
	/bin/rm -f ./out/ueq*.bmp
	/bin/rm -f ./out/u*.bmp
	/bin/rm -f ./out/vor*.bmp
	/bin/rm -f ./out/force*.bmp
	/bin/rm -f ./out/sforce*.bmp
	/bin/rm -f ./out/cuda_rho*.bmp
	/bin/rm -f ./out/cuda_ueq*.bmp
	/bin/rm -f ./out/cuda_u*.bmp
	/bin/rm -f ./out/cuda_vor*.bmp
	/bin/rm -f ./out/cuda_force*.bmp
	/bin/rm -f ./out/cuda_sforce*.bmp
	/bin/rm -f *~ .*.sw*
	/bin/rm -f ./src/*~ ./src/.*.sw*
	/bin/rm -f ./out/slice*.m
	/bin/rm -f ./out/*.raw

clean:
	/bin/rm -f lb3d.exe lb3d a.out slice.exe a.exe* *.stackdump
