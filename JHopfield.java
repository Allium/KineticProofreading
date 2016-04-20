/**
 * @author Daniel Seeto
 * @version 3/27/16
 * 
 */
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.Arrays;
import java.text.DateFormat;
import java.util.Calendar;
public class JHopfield {
	public static void main(String[] args){
		final long startTime = System.nanoTime();
		Calendar cal = Calendar.getInstance();
		DateFormat df = DateFormat.getDateTimeInstance(DateFormat.FULL, DateFormat.MEDIUM);

		double T = 0;
		double U = 0;
		
		double dummyI=0;
		double dummyC=0;

		double HopfieldRate = 0;
		
		double ERR1 = 0;
		double ERR2 = 0;
		
		// Command Line delta input
		double delta1 = Double.parseDouble(args[0]);
		double delta2 = 1;

		// Default hopfield is on
		double hopfield = 1;
		String filekey = "";
		// Input Hopfield on or off
		if (args[1].equalsIgnoreCase("h")) {
			hopfield = 1;
			filekey = "Hopfield";
		}

		else if (args[1].equalsIgnoreCase("n")) {
			hopfield = 0;
			filekey = "Notfield";
		}
		
		
		// ======================================Unprimed Quantities
		int E1 = 1;
		int A1 = 5000;
		int B1 = 0;
		int C1 = 0;
		int E2 = 1;
		int A2 = 5000;
		int B2 = 0;
		int C2 = 0;
		
		int A10 = A1;
		int A20 = A2;

		// ======================================Primed Quantities
		int pA1 = 5000;
		int pB1 = 0;
		int pC1 = 0;
		int pA2 = 5000;
		int pB2 = 0;
		int pC2 = 0;
		
		int pA10 = pA1;
		int pA20 = pA2;

		
		// ======================================Rates unprimed
		double A1B1 = .04;
		double B1A1 = .005;
		double B1C1 = .01;
		double B1A2 = .04 *delta1 * hopfield;
		double C1B1 = B1C1;
		double C1A1 = .04 * hopfield;
		double C1A2 = .005*delta1;
		double A2C1 = .04;
		
		double pA1B1 = A1B1;
		double pB1A1 = C1A2;
		double pB1C1 = B1C1;
		double pB1A2 = C1A1;
		double pC1B1 = C1B1;
		double pC1A1 = B1A2;
		double pC1A2 = B1A1;
		double pA2C1 = A1B1;

		// Timesteps
		double timeMax = 2000 * (A1 + A2 + pA1 + pA2);

		// ======================================Loop Counters
		int t = 0;
		int countA1 = 1;
		int countB1 = 1;
		int countC1 = 1;
		int countA2 = 1;
		
		int countpA1 = 1;
		int countpB1 = 1;
		int countpC1 = 1;
		int countpA2 = 1;
		
		int countE1 = 1;
		int countE2 = 1;


		// ======================================Irreversible counter
		int xB1A2 = 0;
		int xC1A1 = 0;
		
		int xpC1A1 = 0;
		int xpB1A2 = 0;
		
		//==Products
		int xC1A2  =0;
		int xB1A1 = 0;
		
		int xpC1A2 = 0;
		int xpB1A1 = 0;

		// Sums of irreversibles
		int xu = xB1A2 + xC1A1 + xpC1A1 + xpB1A2;

		// ======================================Entropy
		double i = -A1 * Math.log(A1);
		double j = -A2 * Math.log(A2);
		double k = -pA1 * Math.log(pA1);
		double l = -pA2 * Math.log(pA2);

		double S = i + j + k + l;

		// Initial Entropy
		System.out.println("S =" + S);
		double S0 = S;
		// For change in entropy calculations
		double dS = 0;
		double R = S;
		// Total work and change in work
		int dW = 0;
		int W = 0;
		
		//Incorrect, ∆ incorrect, correct, ∆incorrect, 
		double I = 0;
		double dI = 0;
		double C = 0;
		double dC = 0;
		double dIdC = 0;
		
		//Incorrect/ Correct
		double ratio = 0;
		
		double sum = 0;

		// =================================================================Generate
		// datafile with header
		String fileName = "V2" + filekey + "_" + args[0]+ "_.txt";
		try {
			PrintWriter outputStream = new PrintWriter(fileName);
			outputStream.println(df.format(cal.getTime()));
			outputStream.println("∆ = " + Double.toString(delta1));
			outputStream.println();
			outputStream.println("A1B1\tB1A1\tB1C1\tB1A2\tC1A1\tC1B1\tC1A2\tA2C1");
			outputStream.println(A1B1 + "\t" + B1A1 + "\t" + B1C1 + "\t" + B1A2 + "\t" + C1A1 + "\t" +C1B1 + "\t" + C1A2 + "\t" + A2C1);
			outputStream.println();
			outputStream.println("A1'B1'\tB1'A1'\tB1'C1'\tB1'A2'\tC1'A1'\tC1'B1'\tC1'A2'\tA2'C1'");
			outputStream.println(pA1B1 + "\t" + pB1A1 + "\t" + pB1C1 + "\t" + pB1A2 + "\t" + pC1A1 + "\t" + pC1B1 + "\t" + pC1A2 + "\t" + pA2C1);
			outputStream.println();

			// ======================================Show T, S, W, D+D',
			outputStream.println("Time\tA1\tA2\tA1'\tA2'\tS\tW\tProd\tC\tI\tA1xA2\tA2xA1\tA1'xA2'\tA2'xA1'");

			// Loop through Time and output line by line
			for (t = 0; t <= timeMax; t++) {
				if (t % (timeMax / 1000) == 0) {
					outputStream.print(t);
					outputStream.print("\t");
					 outputStream.print(A1);
					 outputStream.print("\t");
					 outputStream.print(A2);
					 outputStream.print("\t");
					 outputStream.print(pA1);
					 outputStream.print("\t");
					 outputStream.print(pA2);
					 outputStream.print("\t");
					
					// ======================================Entropy
					outputStream.print((S - S0));
					outputStream.print("\t");
					
					// ======================================Work
					outputStream.print(-1 * (W));
					outputStream.print("\t");
					
					// ======================================Product or box
					// transitions
					outputStream.print(-1 * (xB1A2 + xC1A2 + xB1A1 + xC1A1 + xpB1A2 + xpC1A2+ xpB1A1 + xpC1A1));
					outputStream.print("\t");
					//==========================================correct and incorrect
		
					outputStream.print(-1 * (xB1A2 + xC1A2 + xpB1A1 + xpC1A1));
					outputStream.print("\t");
					outputStream.print(-1 * (xB1A1 + xC1A1 + xpB1A2 + xpC1A2));
					outputStream.print("\t");
					outputStream.print(xB1A2 + xC1A2);
					outputStream.print("\t");
					outputStream.print(xB1A1 + xC1A1);
					outputStream.print("\t");
					outputStream.print(xpB1A2 + xpC1A2);
					outputStream.print("\t");
					outputStream.print(xpB1A1 + xpC1A1);
					
					
/*
				//Incorrect / Correct
					if (C==0){
						outputStream.print(0);
					}
					else if (C!=0){
						outputStream.print(I/C);
					}
					outputStream.print("\t");
					

				//Change in Incorrect/Correct
			//		outputStream.print((I/C - ratio)*A1/A2);

			//		outputStream.print("\t");
					
				//Change in top and bottom
					outputStream.print(HopfieldRate);

*/
					outputStream.println();
				}

				// Reset counters after each time step
				
				//Reset Rate
				
				A1B1 = .01 * A1 / A10;
				A2C1 = .01 * A2 / A20;
				
				pA1B1 = .01 * pA1 / pA10;
				pA2C1 = .01 * pA2 / pA20;

				countA1 = 1;
				countA2 = 1;
				countB1 = 1;
				countC1 = 1;
				
				countpA1 = 1;
				countpA2 = 1;
				countpB1 = 1;
				countpC1 = 1;

				
				while (countA1 <= 1 && A1 >0 && (B1+C1)<1){
					double chance = (double)(Math.random());
					if (chance < A1B1){
						A1--;
						B1++;
						T = 1;
					}
					countA1++;
				}

				while (countB1 <= B1 && B1 > 0) {
					double chance = (double) (Math.random());
					if (chance < B1A1) {
						B1--;
						A1++;
						if (T == 0)
							xB1A1++;
						
					} else if (chance >= 1 - B1C1) {
						B1--;
						C1++;
					} else if (chance < (.5+B1A2) && chance > .5){
						B1--;
						A2++;
						if (T == 1)
							xB1A2++;
					}
					countB1++;
				}

				while (countC1 <= C1 && C1 > 0) {
					double chance = (double) (Math.random());
					if (chance < C1A1) {
						C1--;
						A1++;
						if (T ==0)
							xC1A1++;
						
					} else if (chance >= (1 - C1B1)) {
						C1--;
						B1++;
					} else if (chance > .5 && chance <= (.5 + C1A2)) {
						C1--;
						A2++;
						if (T==1)
							xC1A2++;
							
					}
					countC1++;
				}
				
				while (countA2 <= 1 & A2 >0 && (B1+C1)<1){
					double chance = (double) (Math.random());
					if (chance < A2C1){
						A2--;
						C1++;
						T = 0;
					}
					countA2++;
				}

	
				
				while (countpA1 <= 1 && pA1 >0 && (pB1+pC1)<1){
					double chance = (double)(Math.random());
					if (chance < pA1B1){
						pA1--;
						pB1++;
						U = 1;
					}
					countpA1++;
				}

				while (countpB1 <= pB1 && pB1 > 0) {
					double chance = (double) (Math.random());
					if (chance < pB1A1) {
						pB1--;
						pA1++;
						if (U == 0)
							xpB1A1++;
						
					} else if (chance >= 1 - pB1C1) {
						pB1--;
						pC1++;
					} else if (chance < (.5+pB1A2) && chance > .5){
						pB1--;
						pA2++;
						if (U == 1)
							xpB1A2++;
					}
					countpB1++;
				}

				while (countpC1 <= pC1 && pC1 > 0) {
					double chance = (double) (Math.random());
					if (chance < pC1A1) {
						pC1--;
						pA1++;
						if (U ==0)
							xpC1A1++;
						
					} else if (chance >= (1 - pC1B1)) {
						pC1--;
						pB1++;
					} else if (chance > .5 && chance <= (.5 + pC1A2)) {
						pC1--;
						pA2++;
						if (U==1)
							xpC1A2++;
							
					}
					countpC1++;
				}
				
				while (countpA2 <= 1 & pA2 >0 && (pB1+pC1)<1){
					double chance = (double) (Math.random());
					if (chance < pA2C1){
						pA2--;
						pC1++;
						U = 0;
					}
					countpA2++;
				}


				// Calculate entropy, checking for condition 0's that might
				// throw off the calculation

				if (A1 > 0) {
					i = -A1 * Math.log(A1);
				} else if (A1 <= 0) {
					i = 0;
				}

				if (A2 > 0) {
					j = -A2 * Math.log(A2);
				} else if (A2 <= 0) {
					j = 0;
				}
				
				if (pA1 > 0) {
					k = -pA1 * Math.log(pA1);
				} else if (pA1 <= 0) {
					k = 0;
				}

				if (pA2 > 0) {
					l = -pA2 * Math.log(pA2);
				} else if (pA2 <= 0) {
					l = 0;
				}	
				

				// Update ∆S, S, ∆W, W, etc.

				dS = i + j + k + l - R;
				S = i + j + k + l;
				R = S;

				// Update number of product
				
				
				dW = xu - W;
				xu = xB1A2 + xC1A1 + xpB1A2 + xpC1A1;
				W = xu;
				
				
				dI = xB1A1 + xC1A1 + xpC1A2 + xpB1A2- I;
				I = xB1A1 + xC1A1 + xpC1A2 + xpB1A2;
				
				dC = xC1A2 + xB1A2 + xpB1A1 + xpC1A2 - C;
				C = xC1A2 + xB1A2 + xpB1A1 + xpC1A2;
				
				if (dC != 0){
					dIdC = dI / dC;
				}
				else if (dC == 0){
					dIdC = 0;
				}


				//Update the incorrect/correct and the difference
				if (t % (timeMax / 1000) == 0 && C !=0){
					dummyI = I;
					dummyC = C;
					ratio = I/C;
				}
				else if (C ==0){
					ratio = 0;
				}

				if ((C-dummyC > 0)){
					HopfieldRate = (I-dummyI)*A1/A2/(C - dummyC);
				}


				if (((xB1A2 + xC1A2)/A1) > 0){ 
					ERR1 = ( (xpC1A2 + xpB1A2) / pA1 ) / ( (xB1A2 + xC1A2) /A1 );
}

				if (((xpC1A1 + xpB1A1)/pA2) > 0){
					ERR2 = ( (xB1A1 + xC1A1) / A2 ) / ( (xpC1A1 + xpB1A1) /pA2);
}

				
			}
			outputStream.close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}

		// ======================================Terminal
		// Output================================================
		System.out.println("S-S(0) = " + (S - S0));
		// System.out.println();
		// ======================================Particle numbers
		System.out.println("A1\tB1\tC1\tA2");
		System.out.println(A1 + "\t" + B1 + "\t" + C1 + "\t" + A2);
		System.out.println();
		System.out.println("Irreversibles:\t" + xu);
		System.out.println();
		System.out.println("Correct and Incorrect");
		System.out.println("xB1A2\txC1A2\txB1A1\txC1A1");
		System.out.println(xB1A2 + "\t" + xC1A2 + "\t" + xB1A1 + "\t" + xC1A1);
		System.out.println();
		final long duration = System.nanoTime() - startTime;
		System.out.println();
		System.out.println("Execution time: " + duration / 1000000000.0 + " seconds.");
		System.out.println("Saved to " + fileName);
		System.out.println("");
	}
}
