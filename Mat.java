public class Mat {
	
	/**
	 * optional class with static matrix helper functions
	 * directly operate on 2d arrays
	 * mostly derived from http://math.nist.gov/javanumerics/jama/ 
	 */
	
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

}
