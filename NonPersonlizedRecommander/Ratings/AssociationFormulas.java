package Ratings;

import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;

import Ratings.FileReadingAndProcessing;

public class AssociationFormulas {

	FileReadingAndProcessing fileReading;

	Map<String, Integer> referenceMovie1_SimpleFormulas = new HashMap<String, Integer>();
	Map<String, Integer> referenceMovie2_SimpleFormulas = new HashMap<String, Integer>();
	Map<String, Integer> referenceMovie3_SimpleFormulas = new HashMap<String, Integer>();

	Map<String, Integer> referenceMovie1_AdvancedFormulas = new HashMap<String, Integer>();
	Map<String, Integer> referenceMovie2_AdvancedFormulas = new HashMap<String, Integer>();
	Map<String, Integer> referenceMovie3_AdvancedFormulas = new HashMap<String, Integer>();

	// to store the temporary values of (X and Y)/(!X and Y)
	Map<String, Double> tempRatioMovie1_AdvancedFormulas = new HashMap<String, Double>();
	Map<String, Double> tempRatioMovie2_AdvancedFormulas = new HashMap<String, Double>();
	Map<String, Double> tempRatioMovie3_AdvancedFormulas = new HashMap<String, Double>();

	ArrayList<Map<String, Integer>> referenceMovies_SimpleFormulas;
	ArrayList<Map<String, Integer>> referenceMovies_AdvancedFormulas;
	// to store the temporary values of (X and Y)/(!X and Y)
	ArrayList<Map<String, Double>> tempRatioMovies_AdvancedFormulas;

	// # of users that has rated or not rated the ref movie
	int[] ratedAndNotRatedMovieX;

	// top 5 rating associations - simple formula
	int[] topAssMov;

	// top 5 rating associations - advanced formula
	double[] topAssMov2;

	// top 5 rating associations value computed by formulas
	double[] topAssMovVals;

	public AssociationFormulas(int topN, int movieAmount) {
		fileReading = new FileReadingAndProcessing();
		topAssMov = new int[2 * topN * movieAmount];
		topAssMov2 = new double[2 * topN * movieAmount];
		topAssMovVals = new double[4 * topN * movieAmount];
		// referenceMovies_SimpleFormulas = new ArrayList<Map<String,
		// Integer>>();
		// referenceMovies_AdvancedFormulas = new ArrayList<Map<String,
		// Integer>>();
		// tempRatioMovies_AdvancedFormulas = new ArrayList<Map<String,
		// Double>>();
	}

	public void scanForStatistics(String movie1, String movie2, String movie3) {

		// to make it flexible
		// for (int i = 0; i < movies.length; i++) {
		// Map<String, Integer> referenceMovie_SimpleFormulas = new
		// HashMap<String, Integer>();
		// referenceMovies_SimpleFormulas.add(referenceMovie_SimpleFormulas);
		// }
		// for (int i = 0; i < movies.length; i++) {
		// Map<String, Integer> referenceMovie_AdvancedFormulas = new
		// HashMap<String, Integer>();
		// referenceMovies_AdvancedFormulas
		// .add(referenceMovie_AdvancedFormulas);
		// }
		// for (int i = 0; i < movies.length; i++) {
		// Map<String, Double> tempRatioMovie_AdvancedFormulas = new
		// HashMap<String, Double>();
		// tempRatioMovies_AdvancedFormulas
		// .add(tempRatioMovie_AdvancedFormulas);
		// }

		// initialize all the reference lists
		fileReading.creatMovieList(referenceMovie1_SimpleFormulas,
				referenceMovie2_SimpleFormulas, referenceMovie3_SimpleFormulas,
				referenceMovie1_AdvancedFormulas,
				referenceMovie2_AdvancedFormulas,
				referenceMovie3_AdvancedFormulas,
				"recsys-data-movie-titles.csv");

		ratedAndNotRatedMovieX = fileReading.constructAssociationList(
				referenceMovie1_SimpleFormulas, referenceMovie2_SimpleFormulas,
				referenceMovie3_SimpleFormulas,
				referenceMovie1_AdvancedFormulas,
				referenceMovie2_AdvancedFormulas,
				referenceMovie3_AdvancedFormulas, "recsys-data-ratings.csv",
				movie1, movie2, movie3);

	}

	public void sortAndSelect(Map<String, Integer> referenceMovies, int index,
			String key) {
		// search the maximum 5 elements in all association rating lists, remove
		// the key itself
		referenceMovies.remove(key);
		// put the <Key, Value> entryset into a Linkedlist
		List<Entry<String, Integer>> list = new LinkedList<Entry<String, Integer>>(
				referenceMovies.entrySet());

		Collections.sort(list, new Comparator<Object>() {
			// sorting by the values from max to min
			public int compare(Object o1, Object o2) {
				return ((Comparable<Integer>) ((Map.Entry<String, Integer>) (o2))
						.getValue())
						.compareTo(((Map.Entry<String, Integer>) (o1))
								.getValue());
			}
		});

		int i = 0;
		for (Iterator<Entry<String, Integer>> it = list.iterator(); it
				.hasNext();) {
			Map.Entry<String, Integer> entry = it.next();
			topAssMov[index + i] = new Integer(entry.getKey());
			topAssMov[index + i + 1] = entry.getValue();

			if (i == 8) {
				break;
			}
			i = i + 2;
		}
	}

	public void sortAndSelect2(Map<String, Double> referenceMovies, int index,
			String key) {
		// search the maximum 5 elements in all association rating lists,
		// firstly remove
		// the key itself
		referenceMovies.remove(key);
		// put the <Key, Value> entryset into a Linkedlist
		List<Entry<String, Double>> list = new LinkedList<Entry<String, Double>>(
				referenceMovies.entrySet());

		Collections.sort(list, new Comparator<Object>() {
			// sorting by the values from max to min
			public int compare(Object o1, Object o2) {
				return ((Comparable<Double>) ((Map.Entry<String, Double>) (o2))
						.getValue())
						.compareTo(((Map.Entry<String, Double>) (o1))
								.getValue());
			}
		});

		// obtain the top5
		int i = 0;
		for (Iterator<Entry<String, Double>> it = list.iterator(); it.hasNext();) {
			Map.Entry<String, Double> entry = it.next();
			topAssMov2[index + i] = new Integer(entry.getKey());
			topAssMov2[index + i + 1] = entry.getValue();

			if (i == 8) {
				break;
			}
			i = i + 2;
		}
	}

	public void computeAndOutput(String movie1, String movie2, String movie3,
			int tops) {

		int index = 0;
		int interval = tops * 2;

		// to define (X and Y)/(!X and Y)

		this.convert(referenceMovie1_SimpleFormulas,
				referenceMovie1_AdvancedFormulas,
				tempRatioMovie1_AdvancedFormulas);
		this.convert(referenceMovie2_SimpleFormulas,
				referenceMovie2_AdvancedFormulas,
				tempRatioMovie2_AdvancedFormulas);
		this.convert(referenceMovie3_SimpleFormulas,
				referenceMovie3_AdvancedFormulas,
				tempRatioMovie3_AdvancedFormulas);

		// to find the top5 maximum (X and Y)
		this.sortAndSelect(referenceMovie1_SimpleFormulas, index, movie1);
		this.sortAndSelect(referenceMovie2_SimpleFormulas, index + interval,
				movie2);
		this.sortAndSelect(referenceMovie3_SimpleFormulas,
				index + interval * 2, movie3);

		// to find the top5 maximum (X and Y)/(!X and Y)
		this.sortAndSelect2(tempRatioMovie1_AdvancedFormulas, index, movie1);
		this.sortAndSelect2(tempRatioMovie2_AdvancedFormulas, index + interval
				* 1, movie2);
		this.sortAndSelect2(tempRatioMovie3_AdvancedFormulas, index + interval
				* 2, movie3);

		// compute using the simple formula
		for (int i = 0; i < 30; i = i + 2) {
			int j = i;
			topAssMovVals[i] = topAssMov[i];
			double originalValue = (double) topAssMov[i + 1]
					/ (double) ratedAndNotRatedMovieX[0 + (int) j / interval];
			BigDecimal bg = new BigDecimal(originalValue);
			double value = bg.setScale(2, BigDecimal.ROUND_HALF_UP)
					.doubleValue();
			topAssMovVals[i + 1] = value;
		}
		// compute using the advanced formula
		for (int i = 0; i < 30; i = i + 2) {
			int j = i;
			topAssMovVals[i + 30] = topAssMov2[i];

			// (X and Y)/(!X and Y) * (!X / X)
			double originalValue = (double) topAssMov2[i + 1]
					* (double) ratedAndNotRatedMovieX[(int) j / interval + 3]
					/ (double) ratedAndNotRatedMovieX[(int) j / interval];

			BigDecimal bg = new BigDecimal(originalValue);
			double value = bg.setScale(2, BigDecimal.ROUND_HALF_UP)
					.doubleValue();

			// Alternative
			// double value = Math.round(originalValue * 100.0) / 100.0;

			topAssMovVals[i + 31] = value;
		}

		fileReading.writeReseults(topAssMovVals, movie1, movie2, movie3);
	}

	// compute (X and Y)/(!X and Y)
	public void convert(Map<String, Integer> numerators,
			Map<String, Integer> denominators, Map<String, Double> results) {

		Set<Entry<String, Integer>> denoSet = denominators.entrySet();
		Iterator<Entry<String, Integer>> iteratorDeno = denoSet.iterator();

		while (iteratorDeno.hasNext()) {
			Map.Entry<String, Integer> entry = iteratorDeno.next();
			String keyD = entry.getKey();
			int denominator = entry.getValue();

			if (denominator != 0) {
				int numerator = numerators.get(keyD);
				double value = (double) numerator / (double) denominator;
				results.put(keyD, value);
			}
		}

	}

	public static void main(String[] args) {
		int showTop = 5;
		// String[] movies = new String[] { "11", "121", "8587" };
		// String[] movies = new String[] { "7443", "641", "9806" };
		String[] movies = new String[] { "568", "581", "607" };

		AssociationFormulas associationFormulas = new AssociationFormulas(
				showTop, movies.length);
		// associationFormulas.scanForStatistics("11", "121", "8587");
		// associationFormulas.computeAndOutput("11", "121", "8587", showTop);

		// associationFormulas.scanForStatistics("7443", "641", "9806");
		// associationFormulas.computeAndOutput("7443", "641", "9806", showTop);

		associationFormulas.scanForStatistics("568", "581", "607");
		associationFormulas.computeAndOutput("568", "581", "607", showTop);

	}

}
