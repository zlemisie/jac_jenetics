/*
 * Java Genetic Algorithm Library (@__identifier__@).
 * Copyright (c) @__year__@ Franz Wilhelmstötter
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 * Author:
 *    Franz Wilhelmstötter (franz.wilhelmstoetter@gmx.at)
 */
package org.jenetics.example;

import org.jenetics.*;
import org.jenetics.engine.Engine;
import org.jenetics.engine.EvolutionStatistics;
import org.jenetics.util.ISeq;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.*;
import java.util.stream.IntStream;

import static java.lang.Math.*;
import static org.jenetics.engine.EvolutionResult.toBestPhenotype;
import static org.jenetics.engine.limit.bySteadyFitness;

public class TSP2 {

	// Problem initialization:
	// Calculating the adjacence matrix of the "city" distances.

	private static int STOPS = 20;
	//private static final double[][] ADJACENCE = matrix(STOPS);
	//private static final double[][] ADJACENCE = matrix2("tsp/17.txt"); //2085
	//private static final double[][] ADJACENCE = matrix2("tsp/26.txt"); //937
	private static final double[][] ADJACENCE = matrix2("tsp/42.txt"); //699
	//private static final double[][] ADJACENCE = matrix2("tsp/48.txt"); //10628


	private static double[][] matrix2(String plik) {
		int stops;
		double[][] matrix = null;
		try{
			//String curDir = System.getProperty("user.dir");
			//System.out.println(curDir);

			File file = new File(plik);
			Scanner in = new Scanner(file);
			stops = in.nextInt();
			STOPS = stops;

			double liczba;
			matrix = new double[stops][stops];


			System.out.println();
			for (int i = 0; i < stops; ++i) {
				for (int j = 0; j < stops; ++j) {
					//System.out.print("i: "+i);
					liczba = in.nextDouble();
					matrix[i][j] = liczba;
					//System.out.print(liczba + " ");
				}
				//System.out.println();

			}

			



		} catch (FileNotFoundException e){
			System.out.println("FNFE");
		}

		return matrix;
	}

	private static double[][] matrix(int stops) {

		final double radius = 10.0;
		double[][] matrix = new double[stops][stops];

		for (int i = 0; i < stops; ++i) {
			for (int j = 0; j < stops; ++j) {
				matrix[i][j] = chord(stops, abs(i - j), radius);
			}
		}
		return matrix;
	}

	private static double chord(int stops, int i, double r) {
		return 2.0*r*abs(sin((PI*i)/stops));
	}

	/*// Calculate the path length of the current genotype.
	public static
	Double dist(final Genotype<EnumGene<Integer>> gt) {
		// Convert the genotype to the traveling path.
		final int[] path = gt.getChromosome().toSeq().stream()
			.mapToInt(EnumGene<Integer>::getAllele)
			.toArray();

		// Calculate the path distance.
		return IntStream.range(0, STOPS)
			.mapToDouble(i ->
				ADJACENCE[path[i]][path[(i + 1)%STOPS]])
			.sum();
	}*/

	// Calculate the path length of the current genotype.
	public static
	Double dist(final Genotype<EnumGene<Integer>> gt) {
		// Convert the genotype to the traveling path.
		return dist2(gt.getChromosome().toSeq());
	}

	public static Double dist2(ISeq<EnumGene<Integer>> genes){
		final int[] path = genes.stream()
				.mapToInt(EnumGene<Integer>::getAllele)
				.toArray();

		return IntStream.range(0, STOPS)
				.mapToDouble(i ->
						ADJACENCE[path[i]][path[(i + 1)%STOPS]])
				.sum();
	}

	//20 stops: 62.57378601609235
	public static void main(String[] args) {
		final Engine<EnumGene<Integer>, Double> engine = Engine
				.builder(
						TSP2::dist,
						PermutationChromosome.ofInteger(STOPS))
				.optimize(Optimize.MINIMUM)
						//.maximalPhenotypeAge(100)
				.populationSize(param0)
				.alterers(
						new PartiallyMatchedCrossover<>(param1),
						new TSPMutator1<>(param2)
				)
				.build();

		// Create evolution statistics consumer.
		//final EvolutionStatistics<Double, ?>
				statistics = EvolutionStatistics.ofNumber();

		//final Phenotype<EnumGene<Integer>, Double>
				best = engine.stream()
						// Truncate the evolution stream after 7 "steady"
						// generations.
						.limit(bySteadyFitness(15))
								// The evolution will stop after maximal 100
								// generations.
						.limit(100)
								// Update the evaluation statistics after
								// each generation
						.peek(statistics)
								// Collect (reduce) the evolution stream to
								// its best phenotype.
						.collect(toBestPhenotype());

		System.out.println(statistics);
		System.out.println(best);

		//System.out.println("best:" + best.getFitness());
	}

	static EvolutionStatistics<Double, ?> statistics;
	static Phenotype<EnumGene<Integer>, Double> best;
	private static int param0 = 500;
	private static double param1 = 0.80;
	private static double param2 = 0.20;


	private class ParamObject{

		public double meanValue;
		public LinkedList<Double> values= new LinkedList<>();

		public int p0;
		public double p1;
		public double p2;

		public ParamObject(int p0, double p1, double p2) {
			this.p0 = p0;
			this.p1 = p1;
			this.p2 = p2;
		}

		public String toString(){
			return "v: " + meanValue + " p0: " + p0 + " p1: " + p1 + " p2: " + p2;
		}

		public void meanValue(){
			meanValue = 0;
			for (double d: values){
				meanValue += d;
			}
			meanValue /= values.size();
		}
	}

	public void parameterTest(int howManyTimesEach){

		LinkedList<ParamObject> poList = new LinkedList<>();

		ArrayList<Integer> lp0 = new ArrayList<Integer>(Arrays.asList(500,1000));
		ArrayList<Double> lp1 = new ArrayList<Double>(Arrays.asList(0.5,0.6,0.7,0.8,0.9,1.0));
		ArrayList<Double> lp2 = new ArrayList<Double>(Arrays.asList(0.0,0.1,0.2,0.3,0.4,0.5));



		for (Integer p0: lp0){
			for(Double p1: lp1){
				for(Double p2: lp2){
					ParamObject po = new ParamObject(p0,p1,p2);
					poList.add(po);
					param0 = po.p0;
					param1 = po.p1;
					param2 = po.p2;

					for(int i=0; i<howManyTimesEach; i++) {

						main(null);

						System.out.println("best:" + best.getFitness());

						po.values.add(best.getFitness());

					}
					po.meanValue();
					System.out.println();
					}
			}

			/*Collections.sort(poList, new Comparator<ParamObject>() {
				@Override
				public int compare(ParamObject po1, ParamObject po2) {
					return Double.compare(po1.value, po2.value);
				}
			});*/


		}

		Collections.sort(poList, (a, b) -> Double.compare(a.meanValue, b.meanValue));

		System.out.println("Wyniki: ");
		for (ParamObject po: poList){
			System.out.println(po.toString());
		}

	}

}
