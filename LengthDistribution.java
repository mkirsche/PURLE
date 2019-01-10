/*
 * Given a distribution of read lengths and genome length, this software outputs the following:
 *   1. The expected distribution of read lengths of the set of reads which are not contained in any other read
 *   2. The expected number of non-contained reads below each of a pre-defined set of length thresholds
 */
import java.util.*;
import java.io.*;
public class LengthDistribution {
	
	/*
	 * The number of reads in the dataset with each length
	 */
	static TreeMap<Integer, Integer> lengthFreq;
	
	/*
	 * For each length, the probability of a read with that length being uncontained
	 */
	static TreeMap<Integer, Double> probs;
	
	static String fn, ofn;
	static int genomeLength;
	
	/*
	 * Outputs the expected number of non-contained reads with length <= each of these thresholds
	 */
	static int[] thresholds = {500, 1000, 2500, 5000, 10000, 20000, 50000};
	
	/*
	 * Whether or not to generate a sample length distribution instead of the expected one
	 */
	static boolean sample = false;
	
public static void main(String[] args) throws IOException
{
	if(args.length != 3 && args.length != 4)
	{
		System.out.println("Usage: java LengthDistribution Inputfile Outputfile Genomelength [--sample]");
		System.out.println("Inputfile should contain the length of one sequence (in bp) per line");
		System.out.println("  Alternately, the input file can be in FASTQ or FASTA format.");
		System.out.println("Outputfile is the file to output sample noncontained read lengths to");
		System.out.println("Genomelength is an integer representing the length (in bp) of the genome");
		System.out.println("--sample is an optional flag to produce a sample length distribution instead of the average");
		return;
	}
	fn = args[0];
	ofn = args[1];
	genomeLength = Integer.parseInt(args[2]);
	if(args.length > 3 && args[3].equals("--sample"))
	{
		sample = true;
	}
	lengthFreq = new TreeMap<Integer, Integer>(Collections.reverseOrder());
	Scanner input = new Scanner(new FileInputStream(new File(fn)));
	boolean fasta = false;
	while(input.hasNext())
	{
		int x = 0;
		String s = input.nextLine();
		if(s.startsWith("@"))
		{
			// Fastq format
			x = input.nextLine().length();
			for(int i = 0; i<2; i++) input.nextLine();
		}
		else if(s.startsWith(">"))
		{
			// Fasta format
			fasta = true;
			continue;
		}
		else if(fasta)
		{
			// Fasta format - already scanned the name of the read
			x = 0;
			while(input.hasNext())
			{
				String line = input.nextLine();
				x += line.length();
				if(line.startsWith(">")) break;
			}
		}
		else
		{
			// List of lengths
			x = input.nextInt();
		}
		if(!lengthFreq.containsKey(x)) lengthFreq.put(x, 0);
		lengthFreq.put(x, lengthFreq.get(x)+1);
	}

	/*
	 * Calculate the probability of each read length being uncontained
	 */
	
	probs = new TreeMap<Integer, Double>();

	for(int x : lengthFreq.keySet())
	{
		double probUncontained = 1;
		for(int y : lengthFreq.keySet())
		{
			if(y == x)
			{
				probUncontained *= Math.pow(1.0 * (genomeLength - 1) / genomeLength, lengthFreq.get(y)-1);
				break;
			}
			else
			{
				probUncontained *= Math.pow(1.0 * (genomeLength - (1 + y - x)) / genomeLength, lengthFreq.get(y));
			}
		}
		probs.put(x, probUncontained);
	}
	
	makeDistribution(sample);
	
	for(int t : thresholds)
	{
		printExpectedNumber(t);
	}
	
}
static void makeDistribution(boolean sample) throws IOException
{
	Random r = new Random();
	TreeMap<Integer, Integer> newList = new TreeMap<Integer, Integer>();
	for(int x : lengthFreq.keySet())
	{
		int y = lengthFreq.get(x);
		double p = probs.get(x);
		for(int i = 0; i<y; i++)
		{
			if((sample && r.nextDouble() < p) || (!sample && i <= p*y - 0.5))
			{
				if(!newList.containsKey(x)) newList.put(x, 1);
				else newList.put(x, newList.get(x) + 1);
			}
		}
	}

	PrintWriter out = new PrintWriter(ofn);
	for(int x : newList.keySet())
	{
		for(int i = 0; i<newList.get(x); i++)
			out.println(x);
	}
		
	out.close();
	
	System.out.println("Output non-contained read lengths to " + ofn);

}
static void printExpectedNumber(int threshold)
{
	long tot = 0;
	double expected = 0;
	for(int x :  lengthFreq.keySet())
	{
		if(x <= threshold)
		{
			tot += lengthFreq.get(x);
			expected += lengthFreq.get(x) * probs.get(x);
		}
	}
	System.out.println();
	System.out.println("Total reads of length <= " + threshold + ": " + tot);
	System.out.println("Expected number of uncontained reads of length <= " + threshold + ": " 
			+ String.format("%.2f", expected));
	System.out.println("Percentage uncontained: " + String.format("%.2f", 100 * expected / tot));
}
}
