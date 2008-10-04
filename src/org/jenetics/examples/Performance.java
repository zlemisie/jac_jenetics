/*
 * Java Genetic Algorithm Library (@!identifier!@).
 * Copyright (c) @!year!@ Franz Wilhelmstötter
 *  
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * 
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Author:
 *     Franz Wilhelmstötter (franz.wilhelmstoetter@gmx.at)
 *     
 */
package org.jenetics.examples;

import static java.lang.Math.cos;
import static java.lang.Math.sin;
import static java.lang.Math.toRadians;

import javax.measure.Measure;
import javax.measure.unit.SI;

import org.jenetics.DoubleChromosome;
import org.jenetics.DoubleGene;
import org.jenetics.FitnessFunction;
import org.jenetics.FitnessScaler;
import org.jenetics.Genotype;
import org.jenetics.GenotypeFactory;
import org.jenetics.IdentityScaler;
import org.jenetics.Phenotype;
import org.jenetics.Population;
import org.jscience.mathematics.number.Float64;

/**
 * @author <a href="mailto:franz.wilhelmstoetter@gmx.at">Franz Wilhelmstötter</a>
 * @version $Id: Performance.java,v 1.1 2008-09-24 20:20:11 fwilhelm Exp $
 */
public class Performance {
	
	private static final class Function implements FitnessFunction<DoubleGene, Float64> {
		private static final long serialVersionUID = 2793605351118238308L;
		
		public Float64 evaluate(final Genotype<DoubleGene> genotype) {
			final DoubleGene gene = genotype.getChromosome().getGene(0);
			final double radians = toRadians(gene.doubleValue());
			return Float64.valueOf(Math.log(sin(radians)*cos(radians)));
		}
	}
	
	public static void main(String[] args) {
		final Function ff = new Function();
		final GenotypeFactory<DoubleGene> gtf = Genotype.valueOf(new DoubleChromosome(0, 360));
		final FitnessScaler<Float64> fs = IdentityScaler.valueOf();
		
		final int size = 1000000;
		final Population<DoubleGene, Float64> population = new Population<DoubleGene, Float64>(size);
		for (int i = 0; i < size; ++i) {
			final Phenotype<DoubleGene, Float64> pt = Phenotype.valueOf(
				gtf.newGenotype(), ff, fs, 0
			);
			population.add(pt);
		}
		
		long start = System.currentTimeMillis();
		for (int i = 0; i < size; ++i) {
			population.get(i).getFitness();
		}
		long stop = System.currentTimeMillis();
		System.out.println(Measure.valueOf(stop - start, SI.MILLI(SI.SECOND)));
		
		start = System.currentTimeMillis();
		for (int i = 0; i < size; ++i) {
			population.get(i).getFitness();
		}
		stop = System.currentTimeMillis();
		System.out.println(Measure.valueOf(stop - start, SI.MILLI(SI.SECOND)));
		
		
	}
	
	
}




