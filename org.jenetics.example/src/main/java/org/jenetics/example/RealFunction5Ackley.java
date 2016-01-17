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

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.LinkedList;

import static org.jenetics.engine.EvolutionResult.toBestPhenotype;
import static org.jenetics.engine.limit.bySteadyFitness;

public class RealFunction5Ackley {

	private static double minV = -5.00;
	private static double maxV = 5.00;

	// This method calculates the fitness for a given genotype.
	public static Double eval(final Genotype<BitGene> gt) {

		final double x = number2(gt.getChromosome(0));
		final double y = number2(gt.getChromosome(1));

		return -20*Math.exp(-0.2*Math.sqrt(0.5*(x*x+y*y)))-Math.exp(0.5*(Math.cos(2*Math.PI*x)+Math.cos(2*Math.PI*y)))+Math.E+20;
	}

	private static int number(Chromosome<BitGene> ch){

		int lastGeneIndex = ch.length()-1;
		int result = 0;

		int multiplier = 1;
		for(int i=lastGeneIndex; i>=0; i--){
			int val = ch.getGene(i).booleanValue()? 1 : 0;
			result +=  val*multiplier;
			multiplier *= 2;
		}

		return result;
	}

	private static double number2(Chromosome<BitGene> ch){

		double r = minV + (double)(number(ch))/(Math.pow(2,ch.length())-1) * (maxV - minV);
		//double r = minV + (double)0/65535 * (maxV - minV);
		return r;
	}

	public static Double eval2(double x, double y){
		return x+y;
	}

	public static Double eval3(DoubleGene x, DoubleGene y){
		return x.doubleValue() + y.doubleValue();
	}

	public static double getMinV() {
		return minV;
	}

	public static double getMaxV() {
		return maxV;
	}

	public static void main(String[] args) {
		//DoubleGene dg = DoubleGene.of(minV, maxV, 1);
		//System.out.println("dg:"+dg.doubleValue());
		//dg = dg.newInstance(-6.0);
		//System.out.println("dg2:"+dg.doubleValue());

		int ChromosomeLength = 24;
		//CharSequence cc = new String("1");
		//BigInteger bi = new BigInteger("256");



		final Engine<BitGene, Double> engine = Engine
			// Create a new builder with the given fitness
			// function and chromosome.
			.builder(
					RealFunction5Ackley::eval,
					//DoubleChromosome.of(0.0, 2.0*PI, 2))
					//DoubleChromosome.of(minV, maxV, 1),
					//DoubleChromosome.of(minV, maxV, 1))
					//BitChromosome.of(ChromosomeLength, 1.00),
					//BitChromosome.of(ChromosomeLength, 0.00))
					BitChromosome.of(ChromosomeLength, 0.50),
					BitChromosome.of(ChromosomeLength, 0.50))
							//BitChromosome.of(5),
							//BitChromosome.of(5))

							.populationSize(param0)
							.optimize(Optimize.MINIMUM)
							.alterers(
									new SinglePointCrossover<>(param1),
									new Mutator_R5A<>(param2)

									//new MeanAlterer<>(0.6)
							)
									// Build an evolution engine with the
									// defined parameters.

							.build();

		// Create evolution statistics consumer.
		//final EvolutionStatistics<Double, ?>
			statistics = EvolutionStatistics.ofNumber();

		//final Phenotype<BitGene, Double>
				best = engine.stream().limit(bySteadyFitness(50))
			//.limit(byExecutionTime(Duration.ofSeconds(8)))
			//.limit(100)

			// Update the evaluation statistics after
			// each generation
			.peek(statistics)
			// Collect (reduce) the evolution stream to
			// its best phenotype.
			.collect(toBestPhenotype());

		//System.out.println(statistics);
		//System.out.println(best);
	}

	static EvolutionStatistics<Double, ?> statistics;
	static Phenotype<BitGene, Double>  best;
	private static int param0 = 500;
	private static double param1 = 0.80;
	private static double param2 = 0.10;


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

		ArrayList<Integer> lp0 = new ArrayList<Integer>(Arrays.asList(100, 200));
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
