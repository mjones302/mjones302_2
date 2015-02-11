matmult version 1.0 1/30/2015

Help can be reached at:
-----------------------
Michael Jones
mjones302@gatech.edu

General Usage Notes
-------------------

-This READ_Me covers the entire program known as matmult. H
-Program that can multiplies matrices in three different ways
 and print out the execution time.
- To change the type of matrix multiplication go into the code and change 
	# define type 1
to either:
	1 - to "naive" multiply
	2 - to use the transpose method
	3 - to use the 2D blocking method with both matrices colmn based
	4 - to use the 2D blocking method with the a matrix transposed to use row-based
- To change the number of threads change the defined variable numThreads

To Compile:
----------
If you are using gcc enter 
----------------------------------------------------
gcc matmult.c -std=c99 -on -lm
----------------------------------------------------
to compile the HomeWork_2 program. 

To Run:
Two methods:
 Method 1:
 -------
 using gcc
 --------- 
 enter : n n1 n2 n3
-This method will cause the code to randomly compile and run the multiply and print out only the
 execution time
 
Method 21:
 -------
 using gcc
 --------- 
 enter : n n1 n2 n3 mA.txt mB.txt mC.txt
-This method will cause the code to multiply the column based matrices mA with mB and print out the
 execution time and the answer into mC.txt

Warnings:
---------
-The transpose and block methods only work for square matrices.
