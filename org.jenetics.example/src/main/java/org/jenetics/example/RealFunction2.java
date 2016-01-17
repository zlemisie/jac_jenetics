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

import static java.lang.Math.*;
import static org.jenetics.engine.EvolutionResult.toBestPhenotype;
import static org.jenetics.engine.limit.bySteadyFitness;

public class RealFunction2 {

	// This method calculates the fitness for a given genotype.
	private static Double eval(final Genotype<DoubleGene> gt) {
		//System.out.println("TEST1: " + gt.getNumberOfGenes());
		//final double x = gt.getGene().doubleValue();
		//final double y = gt.getGene().doubleValue();

		final double x = gt.get(0,0).doubleValue();
		final double y = gt.get(0,1).doubleValue();

		//return cos(0.5 + sin(x))*cos(x);
		return (x-2)*(x-2)+y*y;
	}

	public static void main(String[] args) {
		final Engine<DoubleGene, Double> engine = Engine
			// Create a new builder with the given fitness
			// function and chromosome.
			.builder(
				RealFunction2::eval,
				//DoubleChromosome.of(0.0, 2.0*PI, 2))
				DoubleChromosome.of(-5.0, 5.0, 2))
			.populationSize(1000)
			.optimize(Optimize.MINIMUM)
			.alterers(
				new Mutator<>(0.03),
				new MeanAlterer<>(0.6))
			// Build an evolution engine with the
			// defined parameters.
			.build();

		// Create evolution statistics consumer.
		final EvolutionStatistics<Double, ?>
			statistics = EvolutionStatistics.ofNumber();

		final Phenotype<DoubleGene, Double> best = engine.stream()
			// Truncate the evolution stream after 7 "steady"
			// generations.
			.limit(bySteadyFitness(7))
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
	}
}