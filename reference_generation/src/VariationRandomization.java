/**
* Author: Jan Schroeder
* Modified by Marek Cmero
**/

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Enumeration;
import java.util.HashMap;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.Map;




class RandomVariation extends Variation {


	private static enum SIZE {
		//SMALL,
		MEDIUM, LARGE, GRANDE
	};

	private static Hashtable<SIZE, int[]> size_ranges = new Hashtable<RandomVariation.SIZE, int[]>();
	static {
		//size_ranges.put(SIZE.SMALL, new int[] { 10, 100 });
		size_ranges.put(SIZE.MEDIUM, new int[] { 300, 1000 });
		size_ranges.put(SIZE.LARGE, new int[] { 2000, 10000 });
		size_ranges.put(SIZE.GRANDE, new int[] { 20000, 100000 });
	}
	private static char[] int_2_seq = { 'A', 'C', 'G', 'T' };

	public RandomVariation(GenomicInterval range, ReferenceIndex chromosomeIndex, boolean homozygousMode,  HashMap<Variation.TYPE, Double> weights) {
		super();

		if(homozygousMode){
			isHomozygous = true;
		}
				/* in het mode, alleles always heterozygous
				 * else {
			double coin = Math.random();
			if (coin < 0.5){
		isHomozygous = true;
			} else
		isHomozygous = false;
		}*/

				TYPE[] vars = TYPE.values();
				double totalWeight = 0.0d;
				Iterator<Map.Entry<Variation.TYPE,Double>> it = weights.entrySet().iterator();
				while (it.hasNext())
				{
			Map.Entry<Variation.TYPE,Double> i = it.next();
			totalWeight += i.getValue();
				}

				int randomIndex = -1;
				double random = Math.random() * totalWeight;

				Iterator<Map.Entry<Variation.TYPE,Double>> it2 = weights.entrySet().iterator();
				double current = 0;
				for(int i=0; i <= vars.length - 1; ++i)
				{
					Double w = weights.get(vars[i]);
					current += w;

					if (random < current)
					{
						randomIndex = i;
						break;
					}
				}

		//int type_index = (int) (Math.random() * (TYPE.values().length - 1)); // exclude SNPs
				TYPE type = TYPE.values()[0];
				try {
					type = TYPE.values()[randomIndex];
				}
				catch (ArrayIndexOutOfBoundsException e) {
					System.out.println("Index out of bounds!");
				}

		int size_index = (int) (Math.random() * SIZE.values().length);
		SIZE size_category = SIZE.values()[size_index];
		int[] size_boundaries = size_ranges.get(size_category);
		int size = Math
		.min((int) (Math.random() * (size_boundaries[1] - size_boundaries[0]))
				+ size_boundaries[0], range.end - range.start);

		int interval_length = range.end - range.start - size;
		int start_coordinate = (int) ((Math.random() * interval_length) + range.start);
		int end_coordinate;
		String insert_sequence = null;

		switch (type) {
		case INSERTION:
			//was in size category small else clause
			int chr_index = (int) (Math.random() * chromosomeIndex.size());
			String chr = chromosomeIndex.getChromosomeName(chr_index);

			int insert_start = (int) (Math.random() * (chromosomeIndex.getLength(chr_index) - size));
			insert_sequence = chr + ":" + insert_start + "-"
			+ (insert_start + size - 1);
			end_coordinate = start_coordinate;
			break;
//			if (size_category == SIZE.SMALL) {
//			// randomize sequence
//			char[] seq = new char[size];
//			for (int i = 0; i < seq.length; i++) {
//		int random_char_index = (int) (Math.random() * 4);
//		seq[i] = int_2_seq[random_char_index];
//			}
//			insert_sequence = new String(seq);
//		} else {
//		}
		case TRANSLOCATION:
			int locationIndex = chromosomeIndex.getIndex(range.chrom);
			chr_index = (int) (Math.random() * (chromosomeIndex.size() - locationIndex)) + locationIndex;
			chr = chromosomeIndex.getChromosomeName(chr_index);

			insert_start = (int) (Math.random() * (chromosomeIndex.getLength(chr_index) - size));
			insert_sequence = chr + ":" + insert_start + "-"
			+ (insert_start + size - 1);
			end_coordinate = start_coordinate;
			break;
		case TANDEM:
			chr_index = chromosomeIndex.getIndex(range.chrom);
			if ((start_coordinate + 2 * size - 1) > chromosomeIndex.getLength(chr_index)) {
		size = (chromosomeIndex.getLength(chr_index) - start_coordinate) / 2;
			}
			int tandem_start = start_coordinate + size;
			insert_sequence = range.chrom + ":" + tandem_start + "-"
			+ (tandem_start + size - 1);
			end_coordinate = start_coordinate + size - 1;
			break;
		case SNP:
			chr_index = (int) (Math.random() * chromosomeIndex.size());
			chr = chromosomeIndex.getChromosomeName(chr_index);

			insert_start = (int) (Math.random() * (chromosomeIndex.getLength(chr_index) - size));
			int random_char_index = (int) (Math.random() * 4);
			insert_sequence = "" + int_2_seq[random_char_index];
			end_coordinate = start_coordinate;
			break;
		default:
			end_coordinate = start_coordinate + size - 1;
		}

		this.setType(type);
		this.setLocation(new GenomicInterval(range.chrom, new Interval(
		start_coordinate, end_coordinate)));
		this.setSequence(insert_sequence);
	}

}

public class VariationRandomization {

	/**
	 * @param args
	 * @throws IOException
	 */
	public static void main(String[] args) throws IOException {

		if (args.length < 4) {
			System.out
			.println("Usage: <input reference index> <variation bin size> <output prefix> <mode (het/hom)> <weights>");
			System.exit(0);
		}

		boolean homozygous_mode = true;
		if (args[3].equals("het"))
			homozygous_mode = false;

		FileWriter outfile1, outfile2=null;
		if(homozygous_mode){
			outfile1 = new FileWriter(args[2]+".fa");
		} else {
			outfile1 = new FileWriter(args[2]+"_allele1.fa");
			outfile2 = new FileWriter(args[2]+"_allele2.fa");
		}
		HashMap<Variation.TYPE, Double> weights = new HashMap<Variation.TYPE, Double>();
				weights.put(Variation.TYPE.INSERTION, 0.0);
				weights.put(Variation.TYPE.DELETION, 0.0);
				weights.put(Variation.TYPE.INVERSION, 0.0);
				weights.put(Variation.TYPE.TANDEM, 0.0);
				weights.put(Variation.TYPE.TRANSLOCATION, 0.0);
				weights.put(Variation.TYPE.SNP, 0.0);

				if (args.length > 4) {
						String[] weightsString = args[4].split(",");
						for (int i=0; i<weightsString.length; i++) {
								String[] weight = weightsString[i].split("=");
								String type = weight[0];
								double val   = Double.parseDouble(weight[1]);
								switch(type) {
										case "dup": weights.put(Variation.TYPE.TANDEM, val);
																break;
										case "del": weights.put(Variation.TYPE.DELETION, val);
																break;
										case "inv": weights.put(Variation.TYPE.INVERSION, val);
																break;
										case "trx": weights.put(Variation.TYPE.TRANSLOCATION, val);
																break;
										case "snp": weights.put(Variation.TYPE.SNP, val);
																break;
										case "ins": weights.put(Variation.TYPE.INSERTION, val);
																break;
								}
						}
				} else {
						weights.put(Variation.TYPE.INSERTION, 0.0);
						weights.put(Variation.TYPE.DELETION, 0.4);
						weights.put(Variation.TYPE.INVERSION, 0.3);
						weights.put(Variation.TYPE.TANDEM, 0.2);
						weights.put(Variation.TYPE.TRANSLOCATION, 0.1);
						weights.put(Variation.TYPE.SNP, 0.0);
				}

		int bucket_size = Integer.parseInt(args[1]);
		ReferenceIndex chromosomeIndex = new ReferenceIndex(args[0]);

		Variation[][] all_variations = new Variation[chromosomeIndex.size()][];


		for (int i=0; i< chromosomeIndex.size(); i++) {

			int no_buckets = (int) Math.ceil((double)chromosomeIndex.getLength(i)/bucket_size);
			all_variations[i] = new Variation[no_buckets];
			System.out.println("Creating "+no_buckets+" variations for chromosome "+chromosomeIndex.getChromosomeName(i));
		}



		for (int i = 0; i < chromosomeIndex.size(); i++) {
			for (int j=0; j<all_variations[i].length -1; j++){
		if(all_variations[i][j] != null){
			outfile1.write(all_variations[i][j].toString()+"\n");
			if(!homozygous_mode && all_variations[i][j].isHomozygous){
				outfile2.write(all_variations[i][j].toString()+"\n");
			}
			continue;
		}
		RandomVariation v;

		while(true){
			v = new RandomVariation(
				new GenomicInterval(
				chromosomeIndex.getChromosomeName(i),
				new Interval(j*bucket_size, Math.min((j+1)*bucket_size, chromosomeIndex.getLength(i)))),
				chromosomeIndex, homozygous_mode, weights);
			if(v.getType() != Variation.TYPE.TRANSLOCATION){
				break;
			}
			GenomicInterval translocated_locus = v.parseSequenceAsInterval();
			int chr_index = chromosomeIndex.getIndex(translocated_locus.chrom);
			int bucket_index = translocated_locus.start/bucket_size;
			if(chr_index != i || bucket_index != j){
				all_variations[chr_index][bucket_index] = new Variation(Variation.TYPE.DELETION, null, translocated_locus);
				all_variations[chr_index][bucket_index].setHomozygous(v.isHomozygous);
				break;
			}
		}
		all_variations[i][j] = v;
		outfile1.write(all_variations[i][j].toString()+"\n");
		if(!homozygous_mode && all_variations[i][j].isHomozygous){
			outfile2.write(all_variations[i][j].toString()+"\n");
		}
			}
		}

		outfile1.flush();
		outfile1.close();
		if(! homozygous_mode){
			outfile2.flush();
			outfile2.close();
		}

	}

}
