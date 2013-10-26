public class LUD {
	
	/**
	 * LU Decomposition
	 * derived from http://math.nist.gov/javanumerics/jama/
	 * (removed dependency to the rest of the library)
	 */

	private double[][] LU; // internal array storage
	private int m, n, pivsign; // matrix dimensions and pivot sign
	private int[] piv; // pivot vector

	public LUD(double[][] M) {

		// Use a "left-looking", dot-product, Crout/Doolittle algorithm.
		LU = M;
		m = M.length;
		n = M[0].length;
		piv = new int[m];
		for (int i = 0; i < m; i++) {
			piv[i] = i;
		}
		pivsign = 1;
		double[] LUrowi;
		double[] LUcolj = new double[m];

		// Outer loop.
		for (int j = 0; j < n; j++) {

			// Make a copy of the j-th column to localize references.
			for (int i = 0; i < m; i++) {
				LUcolj[i] = LU[i][j];
			}

			// Apply previous transformations.
			for (int i = 0; i < m; i++) {
				LUrowi = LU[i];

				// Most of the time is spent in the following dot product.
				int kmax = Math.min(i, j);
				double s = 0.0;
				for (int k = 0; k < kmax; k++) {
					s += LUrowi[k] * LUcolj[k];
				}

				LUrowi[j] = LUcolj[i] -= s;
			}

			// Find pivot and exchange if necessary.
			int p = j;
			for (int i = j + 1; i < m; i++) {
				if (Math.abs(LUcolj[i]) > Math.abs(LUcolj[p])) {
					p = i;
				}
			}
			if (p != j) {
				for (int k = 0; k < n; k++) {
					double t = LU[p][k];
					LU[p][k] = LU[j][k];
					LU[j][k] = t;
				}
				int k = piv[p];
				piv[p] = piv[j];
				piv[j] = k;
				pivsign = -pivsign;
			}

			// Compute multipliers.
			if (j < m & LU[j][j] != 0.0) {
				for (int i = j + 1; i < m; i++) {
					LU[i][j] /= LU[j][j];
				}
			}
		}
	}

	public boolean isNonsingular() {
		for (int j = 0; j < n; j++) {
			if (LU[j][j] == 0)
				return false;
		}
		return true;
	}

	public int[] getPivot() {
		int[] p = new int[m];
		for (int i = 0; i < m; i++) {
			p[i] = piv[i];
		}
		return p;
	}   

	// Return pivot permutation vector as a one-dimensional double array
	public double[] getDoublePivot() {
		double[] vals = new double[m];
		for (int i = 0; i < m; i++) {
			vals[i] = (double) piv[i];
		}
		return vals;
	}
	
	public double det() {
		if (m != n) {
			throw new IllegalArgumentException("Matrix must be square.");
		}
		double d = (double) pivsign;
		for (int j = 0; j < n; j++) {
			d *= LU[j][j];
		}
		return d;
	}

	public static String str(double[][] M) {
		StringBuilder sb = new StringBuilder();
		for (int i = 0; i < M.length; i++) {
			for (int j = 0; j < M[i].length; j++) {
				sb.append('\t').append(M[i][j]);
			}
			sb.append('\n');
		}
		return sb.toString();
	}
	
	public static double[][] identity(int m) {
		double[][] I = new double[m][m];
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < m; j++) {
				I[i][j] = i == j ? 1 : 0;
			}
		}
		return I;
	}
	
	
	public static double[][] multiply(double[][] A, double[][] B) {
		int Am = A.length, An = A[0].length;
		int Bm = B.length, Bn = B[0].length;

		if (Bm != An) {
			throw new IllegalArgumentException("Matrix inner dimensions must agree.");
		}
		double[][] C = new double[Am][Bn];
		double[] Bcolj = new double[An];
		for (int j = 0; j < Bn; j++) {
			for (int k = 0; k < An; k++) {
				Bcolj[k] = B[k][j];
			}
			for (int i = 0; i < Am; i++) {
				double[] Arowi = A[i];
				double s = 0;
				for (int k = 0; k < An; k++) {
					s += Arowi[k] * Bcolj[k];
				}
				C[i][j] = s;
			}
		}
		return C;
	}	
	
	public double[][] inverse() {
		return solve(identity(m));
	}
	
	// we need this, as java's clone() method does not perform a deep copy on 2d arrays
	public static double[][] copy(double[][] A) {
		int m = A.length, n = A[0].length; 
		double[][] B = new double[m][n];
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				B[i][j] = A[i][j];
			}
		}
		return B;
	}
   
	public double[][] submat(double[][] A, int[] r, int j0, int j1) {
		double[][] B = new double[r.length][j1 - j0 + 1];
		for (int i = 0; i < r.length; i++) {
			for (int j = j0; j <= j1; j++) {
				B[i][j - j0] = A[r[i]][j];
			}
		}
		return B;
	}
	
	public double[][] getL() {
		double[][] L = new double[m][n];
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				if (i > j) {
					L[i][j] = LU[i][j];
				} else if (i == j) {
					L[i][j] = 1.0;
				} else {
					L[i][j] = 0.0;
				}
			}
		}
		return L;
	}

	public double[][] getU() {
		double[][] U = new double[m][n];
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				if (i <= j) {
					U[i][j] = LU[i][j];
				} else {
					U[i][j] = 0.0;
				}
			}
		}
		return U;
	}	

	/** Solve A*X = B
	@param  B   A Matrix with as many rows as A and any number of columns.
	@return     X so that L*U*X = B(piv,:)
	@exception  IllegalArgumentException Matrix row dimensions must agree.
	@exception  RuntimeException  Matrix is singular.
	*/
	public double[][] solve(double[][] B) {
		if (B.length != m) {
			throw new IllegalArgumentException("Matrix row dimensions must agree.");
		}
		if (!this.isNonsingular()) {
			throw new RuntimeException("Matrix is singular.");
		}

		// Copy right hand side with pivoting
		int nx = B[0].length;
		double[][] X = submat(B, piv, 0, nx - 1);

		// Solve L*Y = B(piv,:)
		for (int k = 0; k < n; k++) {
			for (int i = k + 1; i < n; i++) {
				for (int j = 0; j < nx; j++) {
					X[i][j] -= X[k][j] * LU[i][k];
				}
			}
		}
		// Solve U*X = Y;
		for (int k = n - 1; k >= 0; k--) {
			for (int j = 0; j < nx; j++) {
				X[k][j] /= LU[k][k];
			}
			for (int i = 0; i < k; i++) {
				for (int j = 0; j < nx; j++) {
					X[i][j] -= X[k][j] * LU[i][k];
				}
			}
		}
		return X;
	}
}