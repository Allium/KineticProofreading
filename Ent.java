/**
 *
 * @author Daniel Seeto
 * @version 1/6/16
 */
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.Arrays;
import java.text.DateFormat;
import java.util.Calendar;

public class Ent {
	public static void main(String[] args) {
		final long startTime = System.nanoTime();
		Calendar cal = Calendar.getInstance();
		DateFormat df = DateFormat.getDateTimeInstance(DateFormat.FULL, DateFormat.MEDIUM);
		
		//Command Line delta input
		double delta1 = Double.parseDouble(args[0]);
		double delta2 = 1/delta1;
		
		//Default hopfield is on
		double hopfield = 1;
		String filekey = "";
		//Input Hopfield on or off
		if (args[1].equalsIgnoreCase("h")) {
			hopfield = 1;
			filekey = "Hopfield";
		}

		else if (args[1].equalsIgnoreCase("n")) {
			hopfield = 0;
			filekey = "Notfield";
		}
		//======================================Unprimed Quantities
		int E1 = 1;
		int A1 = 500;
		int B1 = 0;
		int C1 = 0;
		int E2 = 1;
		int A2 = 500;
		int B2 = 0;
		int C2 = 0;
		
		//======================================Primed Quantities
		int pA1 = 500;
		int pB1 = 0;
		int pC1 = 0;
		int pA2 = 500;
		int pB2 = 0;
		int pC2 = 0;
		
		//======================================Rates unprimed
		double A1B1 = .01;
			double B1A1 = .01 * delta2;
		double B1C1 = .01;
		double C1B1 = .01;
		double C1A2 = .01;
	
		double A2B2 = .01;
		double B2A2 = .01 * delta1;
		double B2C2 = .01;
		double C2B2 = .01;
		double C2A1 = .01;
			double C1A1 = .01 * hopfield * delta2;
		double C2A2 = .01 * hopfield * delta1;
		
		//======================================Rates primed
		double pA1B1 = A1B1;
		double pB1A1 = .01 * delta1;
		double pB1C1 = B1C1;
		double pC1B1 = C1B1;
		double pC1A2 = C1A2;
		double pA2B2 = A2B2;
			double pB2A2 = .01 * delta2;
		double pB2C2 = B2C2;
		double pC2B2 = C2B2;
		double pC2A1 = C2A1;
		double pC1A1 = .01 * hopfield * delta1;
			double pC2A2 = .01 * hopfield * delta2;
		
		//======================================Count new products
		int newE1 = 0;
		int newE2 = 0;
		int newA1 = 0;
		int newA2 = 0;
		int newB1 = 0;
		int newB2 = 0;
		int newC1 = 0;
		int newC2 = 0;
		
		int newpA1 = 0;
		int newpA2 = 0;
		int newpB1 = 0;
		int newpB2 = 0;
		int newpC1 = 0;
		int newpC2 = 0;
		
		
		//Timesteps 
		double timeMax =600000;
		
		// ======================================Loop Counters
		int t = 0;
		int countA1 = 1;
		int countB1 = 1;
		int countC1 = 1;
		int countE1 = 1;
		int countA2 = 1;
		int countB2 = 1;
		int countC2 = 1;
		int countE2 = 1;
		
		int countpA1 = 1;
		int countpB1 = 1;
		int countpC1 = 1;
		int countpA2 = 1;
		int countpB2 = 1;
		int countpC2 = 1;
		
		//======================================Irreversible counter
		int xC1A1 = 0;
		int xC1A2 = 0;
		int xC2A2 = 0;
		int xC2A1 = 0;
		
		int xpC1A1 = 0;
		int xpC1A2 = 0;
		int xpC2A2 = 0;
		int xpC2A1 = 0;
		
		//Sums of irreversibles
		int xu = xC1A1 + xC1A2 + xC2A2 + xC2A1;
		int pxu = xpC1A1 + xpC1A2 + xpC2A2 + xpC2A1;
		
		//======================================Entropy 
		double i = -A1 * Math.log(A1);
		double j = -pA1 * Math.log(pA1);
		double k = -A2 * Math.log(A2);
		double l = -pA2 * Math.log(pA2);
		
		double S = i + j + k + l;
		
		//Initial Entropy
		System.out.println("S =" + S);
		double S0 = S;
		//For change in entropy calculations
		double dS = 0;
		double R = S;
		//Total work and change in work
		int dW = 0;
		int W = 0;

		//=================================================================Generate datafile with header
		String fileName = "Results/Sorting/" + filekey + "_" + args[0] + "_.txt";
		try {
			PrintWriter outputStream = new PrintWriter(fileName);
			outputStream.println(df.format(cal.getTime()));
			//================================================================Currently not showing all data
			//outputStream.println("[A1\tB1\tC1\tA2\tB2\tC2\tA1'\tB1'\tC1'\tA2'\tB2'\tC2']");
			//outputStream.println("["+A1 + "\t"+B1 + "\t" + C1 + "\t" + A2 + "\t"+B2 + "\t" + C2 + "\t" 
				//	+pA1 + "\t"+pB1 + "\t" + pC1 + "\t" + pA2 + "\t"+pB2 + "\t" + pC2+ "]");
			//outputStream.println("[∆1\t∆2] ");
			outputStream.println("∆ =" + Double.toString(delta1));
			outputStream.println();
			//outputStream.println("Time\tA1\tA2\tA1'\tA2'\tStep\tStep'\tS\t\t∆S");
			
			//======================================Show T, S, W, D+D', 
			outputStream.println("Time\tS\t\tW\tProducts");
			
			//Loop through Time and output line by line
			for (t = 0; t <= timeMax; t++) {
				
				if ( t%(timeMax/1000)==0 ){
					outputStream.print(t);
					outputStream.print("\t");
	//				outputStream.print(A1);
	//				outputStream.print("\t");
	//				outputStream.print(A2);
	//				outputStream.print("\t");
	//				outputStream.print(pA1);
	//				outputStream.print("\t");
	//				outputStream.print(pA2);
	//				outputStream.print("\t");
	//				outputStream.print(xu);
	//				outputStream.print("\t");
	//				outputStream.print(pxu);
	//				outputStream.print("\t");
					//outputStream.print(S);
					//outputStream.print("\t");
					//======================================Entropy
					outputStream.print(S-S0);
					outputStream.print("\t\t");
	//				outputStream.print(dS);
	//				outputStream.print("\t\t");
	//				outputStream.print(10*dW);
					//======================================Work
					outputStream.print(-1*(W));
					outputStream.print("\t\t");
					//======================================Product or box transitions
					outputStream.print(-1*(xC1A2+xC2A1+ xpC1A2+xpC2A1));
	//				outputStream.print("\t");
	//				outputStream.print(A1);
	//				outputStream.print("\t");
	//				outputStream.print(A2);
	//				outputStream.print("\t");
	//				outputStream.print(pA1);
	//				outputStream.print("\t");
	//				outputStream.print(pA2);
					outputStream.println();
				}
					
				//Reset counters after each time step
				
				countB1 = 1;
				countC1 = 1;
				countE1 = 1;
				countB2 = 1;
				countC2 = 1;
				countE2 = 1;
				
				countpB1 = 1;
				countpC1 = 1;
				countpB2 = 1;
				countpC2 = 1;
				
				newE1 = 0;
				newE2 = 0;
				newA1 = 0;
				newA2 = 0;
				newB1 = 0;
				newB2 = 0;
				newC1 = 0;
				newC2 = 0;
				
				newpA1 = 0;
				newpA2 = 0;
				newpB1 = 0;
				newpB2 = 0;
				newpC1 = 0;
				newpC2 = 0;
				
				//======================================See what happens to E1 if it exists
				
				while (countE1 <= E1 && (A1 + pA1) > 0 && E1 > 0){
					double chance = (double) (Math.random());
					if (chance < A1B1 && A1 > 0){
						E1--;
						A1--;
						B1++;
					}
					else if (chance >= (1 - pA1B1) && pA1 > 0){
						E1--;
						pA1--;
						pB1++;
					}
					countE1++;
				}
				
				//======================================See what happens to B1 if it exists
				
				while (countB1 <= B1 && B1 > 0){
					double chance = (double) (Math.random());
					if (chance < B1A1){
						B1--;
						A1++;
						E1++;
					}
					else if (chance >= 1 - B1C1){
						B1--;
						C1++;
					}
					countB1++;
				}
				
				//====================================== C1 options
				
				while (countC1 <= C1 && C1 >0){
					double chance = (double) (Math.random());
					if (chance < C1A1){
						C1--;
						A1++;
						E1++;
						xC1A1++;
					}
					else if (chance >= (1- C1B1)){
						C1--;
						B1++;
					}
					else if (chance >.5 && chance <= (.5+C1A2)){
						C1--;
						A2++;
						E1++;
						xC1A2++;
					}
					countC1++;
				}
				
				
				//======================================Options E2
				while (countE2 <= E2 && (A2 + pA2) > 0 && E2 > 0 ){
					double chance = (double) (Math.random());
					if (chance < A2B2 && A2 > 0){
						E2--;
						A2--;
						B2++;
					}
					else if (chance >= (1 - pA2B2) && pA2 > 0){
						E2--;
						pA2--;
						pB2++;
					}
					countE2++;
				}
				
				//Options B2
				
				while (countB2 <= B2 && B2 > 0){
					double chance = (double) (Math.random());
					if (chance < B2A2){
						B2--;
						A2++;
						E2++;
					}
					else if (chance >= (1 - B2C2)){
						B2--;
						C2++;
					}
					countB2++;
				}
				
				//======================================Option C2
				while (countC2 <= C2 && C2 >0){
					double chance = (double) (Math.random());
					if (chance < C2A2){
						C2--;
						A2++;
						E2++;
						xC2A2++;
					}
					else if (chance >= (1- C2B2)){
						C2--;
						B2++;
					}
					else if (chance >.5 && chance <= (.5+C2A1)){
						C2--;
						A1++;
						E2++;
						xC2A1++;
					}
					countC2++;
				}
				
				
				//======================================Option B1'
				while (countpB1 <= pB1 && pB1 > 0){
					double chance = (double) (Math.random());
					if (chance < pB1A1){
						pB1--;
						pA1++;
						E1++;
					}
					else if (chance >= 1 - pB1C1){
						pB1--;
						pC1++;
					}
					countpB1++;
				}
				
				
				//======================================Option C1'
				while (countpC1 <= pC1 && pC1 >0){
					double chance = (double) (Math.random());
					if (chance < pC1A1){
						pC1--;
						pA1++;
						E1++;
						xpC1A1++;
					}
					else if (chance >= (1- pC1B1)){
						pC1--;
						pB1++;
					}
					else if (chance >.5 && chance <= (.5+pC1A2)){
						pC1--;
						pA2++;
						E1++;
						xpC1A2++;
					}
					countpC1++;
				}
				
				
				//======================================Option B2'
				while (countpB2 <= pB2 && pB2 > 0){
					double chance = (double) (Math.random());
					if (chance < pB2A2){
						pB2--;
						pA2++;
						E2++;
					}
					else if (chance >= 1 - pB2C2){
						pB2--;
						pC2++;
					}
					countpB2++;
				}
				
				
			//	======================================Option C2'
				while (countpC2 <= pC2 && pC2 >0){
					double chance = (double) (Math.random());
					if (chance < pC2A2){
						pC2--;
						pA2++;
						E2++;
						xpC2A2++;
					}
					else if (chance >= 1- pC2B2){
						pC2--;
						pB2++;
					}
					else if (chance >.5 && chance <= (.5+pC2A1)){
						pC2--;
						pA1++;
						E2++;
						xpC2A1++;
					}
					countpC2++;
				}
				
				
				//Update number of product
				xu = xC1A1 + xC1A2 + xC2A2 + xC2A1;
				pxu = xpC1A1 + xpC1A2 + xpC2A2 + xpC2A1;
				
				
				//Calculate entropy, checking for condition 0's that might throw off the calculation
				
				if (A1 > 0){
					i = -A1 * Math.log(A1);
					}
				else if (A1 <= 0){
					i = 0;
				}
				
				if (pA1 > 0){
					j = -pA1 * Math.log(pA1);
				}
				else if (pA1 <=0){
					j = 0;
				}
				if (A2 > 0){
					k = -A2 * Math.log(A2);
				}
				else if (A2 <= 0){
					k = 0;
				}
				if (pA2 > 0){
					l = -pA2 * Math.log(pA2);
				}
				else if (pA2 <=0){
					l = 0;
				}
				
				//Update ∆S, S, ∆W, W, etc.
				
				
				dS = i + j + k + l - R;
				S = i + j + k + l;
				R = S;
				
				dW = xu + pxu - W;
				xu = xC1A1 + xC1A2 + xC2A2 + xC2A1;
				pxu = xpC1A1 + xpC1A2 + xpC2A2 + xpC2A1;
				W = xu + pxu;
				
				
			}
			outputStream.close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}

		// ======================================Terminal Output================================================
		System.out.println("S-S(0) = "+ (S- S0));
		//System.out.println();
		//======================================Particle numbers
		System.out.println("A1\tB1\tC1\tE1\tA2\tB2\tC2\tE2\tA1'\tB1'\tC1'\tA2'\tB2'\tC2'");
		System.out.println(A1 + "\t" + B1 + "\t" + C1 +"\t"+E1 +"\t"+ A2 + "\t" + B2 + "\t" + C2 +"\t"+E2 + "\t"+pA1 + "\t" + pB1 + "\t" + pC1 +"\t"+ pA2 + "\t" + pB2 + "\t" + pC2);
		System.out.println();
		System.out.println("Total no prime Box 1:\t" +(A1 + B1 + C1));
		System.out.println("Total no prime Box 2:\t" +(A2 + B2 + C2));
		System.out.println("Total Prime Box 1:\t" + (pA1 + pB1 + pC1));
		System.out.println("Total Prime Box 2:\t" + (pA2 + pB2 + pC2));
		System.out.println();
		System.out.println("Irreversibles:\t"+xu);
		System.out.println("Irreversibles':\t"+pxu);
		System.out.println();
		System.out.println("xC1A1\txC1A2\txC2A2\txC2A1\tpxC1A1\tpxC1A2\tpxC2A2\tpxC2A1");
		System.out.println(xC1A1 + "\t" + xC1A2 + "\t"+ xC2A2 + "\t"+ xC2A1 + "\t" + 
				xpC1A1 + "\t" + xpC1A2 + "\t"+ xpC2A2 + "\t"+ xpC2A1);
		final long duration = System.nanoTime() - startTime;
		System.out.println();
		System.out.println("Execution time: " + duration / 1000000000.0 + " seconds.");
		System.out.println("Saved to " + fileName);
		System.out.println("");
	}
}