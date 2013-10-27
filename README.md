LUD
===

single class to perform LU decomposition on 2d arrays (e.g. allows to invert a matrix or solve a system of linear equations).

i generally recommend to use JAMA (http://math.nist.gov/javanumerics/jama/) for matrix operations in java. however if you just need a minimal implementation of the LU decomposition without beeing depended on the JAMA package, this class might be interesting for you. for example you could just copy the LUD class to processing (http://processing.org/) and use it there right away.

Usage
=====

```java
// solve system of linear equations Ax = b
double[][] A = { 
		{ 1, 2, 3 }, 
		{ 4, 5, 6 }, 
		{ 7, 8, 10 } 
};
System.out.println("A:" + Mat.str(A));

double[][] b = { 
		{ 11 }, 
		{ 12 }, 
		{ 13 } 
};
System.out.println("b:" + Mat.str(b));

double[][] x = new LUD(Mat.copy(A)).solve(b); // "A" should not be modified here!
System.out.println("x:" + Mat.str(x));

// test solution
double[][] b_c = Mat.multiply(A, x);
System.out.println("b_c:" + Mat.str(b_c));

// invert matrix A
double[][] Ainv = new LUD(A).inverse(); // allow A to be modified as we wont need it anymore
System.out.println("A^-1:" + Mat.str(Ainv));
```
