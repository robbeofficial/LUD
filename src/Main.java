public class Main {
	
	public static void main(String[] args) {

		// solve system of linear equations Ax = b
		double[][] A = { 
				{ 1, 2, 3 }, 
				{ 4, 5, 6 }, 
				{ 7, 8, 10 } 
		};
		System.out.println("A:" + LUD.str(A));
		
		double[][] b = { 
				{ 11 }, 
				{ 12 }, 
				{ 13 } 
		};
		System.out.println("b:" + LUD.str(b));
		
		double[][] x = new LUD(LUD.copy(A)).solve(b); // A should not be modified here!
		System.out.println("x:" + LUD.str(x));
		
		// invert matrix A
		double[][] Ainv = new LUD(A).inverse(); // allow A to be modified as we wont need it anymore
		System.out.println("A^-1:" + LUD.str(Ainv));
	}
}
