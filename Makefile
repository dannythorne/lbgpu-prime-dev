# i386 and x86-64 Options
#  -mcpu=cpu-type 
#  -march=cpu-type
#  -mfp- math=unit 
#  -masm=dialect 
#  -mno-fancy-math-387
#  -mno-fp-ret-in-387
#  -msoft-float 
#  -msvr3-shlib
#  -mno-wide-multiply 
#  -mrtd 
#  -malign-double
#  -mpreferred-stack-boundary=num
#  -mmmx 
#  -msse 
#  -msse2
#  -msse3
#  -m3dnow
#  -mthreads 
#  -mno-align-stringops 
#  -minline-all-stringops
#  -mpush-args 
#  -maccumulate-outgoing-args 
#  -m128bit-long-double
#  -m96bit-long-double 
#  -mregparm=num 
#  -momit-leaf-frame-pointer
#  -mms-bitfields
#  -mno-red-zone
#  -mcmodel=code-model
#  -m32 
#  -m64


all:
	gcc -x c -g -w -o lb3d -D__device__='' -D__constant__='' ./src/lbgpu_prime.cu -lm

cuda:
	nvcc -v -o lb3d ./src/lbgpu_prime.cu -lm -arch=sm_13

mpi:
	mpicc -w -DPARALLEL=1 -o lb3d ./src/lb3d_prime.c -lm

warn:
	gcc -pedantic -Wall -Wunused-variable -o lb3d ./src/lb3d_prime.c

nocygwin:
	gcc -mno-cygwin -o lb3d ./src/lb3d_prime.c

slice: new_slice.c lbio.c
	gcc -o slice ./src/new_slice.c

spy: spy.c
	gcc -o spy ./src/spy.c

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
	/bin/rm -f *~ .*.sw*
	/bin/rm -f ./src/*~ ./src/.*.sw*
	/bin/rm -f ./out/slice*.m
	/bin/rm -f ./out/new_slice*.m
	/bin/rm -f ./out/*.raw

clean:
	/bin/rm -f lb3d.exe lb3d a.out slice.exe a.exe* *.stackdump dump*.*
