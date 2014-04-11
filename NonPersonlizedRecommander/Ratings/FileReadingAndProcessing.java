package Ratings;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.Map;

public class FileReadingAndProcessing {

	/*
	 * initialize the Maps of movies by recording the <Key, Value> pairs from
	 * the movie file. Key is the movieId, Value is # of associations.
	 */
	public void creatMovieList(Map<String, Integer> movies1,
			Map<String, Integer> movies2, Map<String, Integer> movies3,
			Map<String, Integer> movies4, Map<String, Integer> movies5,
			Map<String, Integer> movies6, String moviesFileName) {

		File file = new File(moviesFileName);
		BufferedReader reader = null;
		try {
			reader = new BufferedReader(new FileReader(file));
			String tempString = null;
			String[] tempStringArray = null;
			// read file by line until reach the end
			while ((tempString = reader.readLine()) != null) {
				// obtain the string segments : <MovieId, MovieName>
				tempStringArray = tempString.split(",");
				// add elements into the Maps
				movies1.put(tempStringArray[0], 0);
				movies2.put(tempStringArray[0], 0);
				movies3.put(tempStringArray[0], 0);
				movies4.put(tempStringArray[0], 0);
				movies5.put(tempStringArray[0], 0);
				movies6.put(tempStringArray[0], 0);

			}

			reader.close();
		} catch (IOException e) {
			e.printStackTrace();
		} finally {
			if (reader != null) {
				try {
					reader.close();
				} catch (IOException e1) {
				}
			}
		}
	}

	/*
	 * scan all the ratings and record the statistical associations for both
	 * computing formulas
	 */
	public int[] constructAssociationList(Map<String, Integer> movies1,
			Map<String, Integer> movies2, Map<String, Integer> movies3,
			Map<String, Integer> movies4, Map<String, Integer> movies5,
			Map<String, Integer> movies6, String moviesFileName, String key1,
			String key2, String key3) {

		// # of users that has rated or not rated the ref movie
		int[] ratedAndNotRated = new int[6];

		File file = new File(moviesFileName);
		BufferedReader reader = null;
		try {
			reader = new BufferedReader(new FileReader(file));
			String tempString = null;
			String[] tempStringArray = null;

			// to justify whether the data of one user is useful. (e.g. YES if
			// he rates Movie X)
			boolean isUseful1 = false;
			boolean isUseful2 = false;
			boolean isUseful3 = false;

			// UserId in the previous line
			String user = "";
			// the IDs of movies which has been rated by certain user
			ArrayList<String> tempRatings = new ArrayList<String>();

			// read file by line until reach the end
			while ((tempString = reader.readLine()) != null) {
				// obtain the string segments : <UserId, MovieId, Rating>
				tempStringArray = tempString.split(",");
				// do statistics according to the formulas
				if (!user.equals("") && !user.equals(tempStringArray[0])) {
					rearrangment(isUseful1, isUseful2, isUseful3, movies1,
							movies2, movies3, movies4, movies5, movies6,
							ratedAndNotRated, tempRatings);
					isUseful1 = false;
					isUseful2 = false;
					isUseful3 = false;
					tempRatings.clear();
				}

				user = tempStringArray[0];
				String movieId = tempStringArray[1];
				tempRatings.add(tempStringArray[1]);

				if (movieId.equals(key1)) {
					isUseful1 = true;
					ratedAndNotRated[0]++;
				}
				if (movieId.equals(key2)) {
					isUseful2 = true;
					ratedAndNotRated[1]++;
				}
				if (movieId.equals(key3)) {
					isUseful3 = true;
					ratedAndNotRated[2]++;
				}

			}

			// record for the last user
			rearrangment(isUseful1, isUseful2, isUseful3, movies1, movies2,
					movies3, movies4, movies5, movies6, ratedAndNotRated,
					tempRatings);
			isUseful1 = false;
			isUseful2 = false;
			isUseful3 = false;
			tempRatings.clear();

			reader.close();
		} catch (IOException e) {
			e.printStackTrace();
		} finally {
			if (reader != null) {
				try {
					reader.close();
				} catch (IOException e1) {
				}
			}
		}

		return ratedAndNotRated;
	}

	// record all effective ratings
	public void reocrdValues(Map<String, Integer> AssociationMovies,
			ArrayList<String> tempRatings) {
		Iterator<String> iterator = tempRatings.iterator();
		while (iterator.hasNext()) {
			String key = iterator.next();
			int newValue = AssociationMovies.get(key) + 1;
			AssociationMovies.put(key, newValue);
		}
	}

	// when refer to a new user, justify whether the overall
	// ratings given by the former user is needed
	public void rearrangment(boolean isUseful1, boolean isUseful2,
			boolean isUseful3, Map<String, Integer> movies1,
			Map<String, Integer> movies2, Map<String, Integer> movies3,
			Map<String, Integer> movies4, Map<String, Integer> movies5,
			Map<String, Integer> movies6, int[] ratedAndNotRated,
			ArrayList<String> tempRatings) {
		if (isUseful1) {
			reocrdValues(movies1, tempRatings);
		} else {
			reocrdValues(movies4, tempRatings);
			ratedAndNotRated[3]++;
		}
		if (isUseful2) {
			reocrdValues(movies2, tempRatings);
		} else {
			reocrdValues(movies5, tempRatings);
			ratedAndNotRated[4]++;
		}
		if (isUseful3) {
			reocrdValues(movies3, tempRatings);
		} else {
			reocrdValues(movies6, tempRatings);
			ratedAndNotRated[5]++;
		}
	}

	public void writeReseults(double[] result, String key1, String key2,
			String key3) {
		// Write the top 5 movies, one per line, to a text file.
		try {
			PrintWriter writer = new PrintWriter("pa1-result.txt", "UTF-8");

			// adjust the output format
			for (int i = 0; i < result.length; i++) {
				if (i == 0 || i == 30) {
					if (i == 30) {
						writer.println();
						writer.println();
					}
					writer.print(key1 + ",");
				} else if (i == 10 || i == 40) {
					writer.println();
					writer.print(key2 + ",");
				} else if (i == 20 || i == 50) {
					writer.println();
					writer.print(key3 + ",");
				}

				if (i % 2 == 0) {
					// movieId
					writer.print((int) result[i] + ",");
				} else if ((i + 1) % 10 == 0) {
					// association
					// creates a new instance of decimal format and tell it that
					// you want the formatted
					DecimalFormat df = new DecimalFormat("0.00");
					writer.print(df.format(result[i]));
				} else {
					DecimalFormat df = new DecimalFormat("0.00");
					writer.print(df.format(result[i]) + ",");
				}
			}

			writer.close();

		} catch (Exception e) {
			System.out.println(e.getMessage());
		}
	}
}
