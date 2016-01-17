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

import static org.jenetics.engine.EvolutionResult.toBestPhenotype;

public class RealFunctionCopy2 {

	private static double minV = -5.00;
	private static double maxV = 5.00;

	// This method calculates the fitness for a given genotype.
	public static Double eval(final Genotype<BitGene> gt) {
		//System.out.println("TEST1: " + gt.getNumberOfGenes());
		//final double x = gt.getGene().doubleValue();
		//final double y = gt.getGene().doubleValue();

		//final double x = gt.get(0,0).doubleValue();
		//final double y = gt.get(0,1).doubleValue();

		//final double x = gt.get(0,0).doubleValue();
		//final double y = gt.get(1,0).doubleValue();


		final double x = 0.0;
		final double y = 0.0;

		//final DoubleGene x = gt.get(0,0);
		//final DoubleGene y = gt.get(0,1);

		//return cos(0.5 + sin(x))*cos(x);
		//return (x-2)*(x-2)+y*y;

		//return x+y;

		/*System.out.println("TTT");
		System.out.println(gt.get(0,0));
		System.out.println(gt.get(0,1));
		System.out.println(gt.get(0,2));
		System.out.println(gt.get(0,15).booleanValue());

		System.out.println(gt.getChromosome(0).length());
		System.out.println("TTT2");
		System.out.println(number(gt.getChromosome(0)));

		System.out.println(number2(gt.getChromosome(0)));*/

		//return eval2(x,y);

		//return 0.0;

		return number2(gt.getChromosome(0))+number2(gt.getChromosome(1));

		//return eval3(x,y);
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
					RealFunctionCopy2::eval,
					//DoubleChromosome.of(0.0, 2.0*PI, 2))
					//DoubleChromosome.of(minV, maxV, 1),
					//DoubleChromosome.of(minV, maxV, 1))
					//BitChromosome.of(ChromosomeLength, 1.00),
					//BitChromosome.of(ChromosomeLength, 0.00))
					BitChromosome.of(ChromosomeLength, 0.50),
					BitChromosome.of(ChromosomeLength, 0.50))
							//BitChromosome.of(5),
							//BitChromosome.of(5))

							.populationSize(100)
							.optimize(Optimize.MINIMUM)
							.alterers(
									new SinglePointCrossover<>(0.80),
									new Mutator_M7<>(0.10)

									//new MeanAlterer<>(0.6)
							)
									// Build an evolution engine with the
									// defined parameters.

							.build();

		// Create evolution statistics consumer.
		final EvolutionStatistics<Double, ?>
			statistics = EvolutionStatistics.ofNumber();

		final Phenotype<BitGene, Double> best = engine.stream()
			//.limit(bySteadyFitness(20))
			//.limit(byExecutionTime(Duration.ofSeconds(8)))
			//.limit(100)

			// Update the evaluation statistics after
			// each generation
			.peek(statistics)
			// Collect (reduce) the evolution stream to
			// its best phenotype.
			.collect(toBestPhenotype());

		System.out.println(statistics);
		System.out.println(best);
	}
}
