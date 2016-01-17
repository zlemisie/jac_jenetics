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
package org.jenetics;

import org.jenetics.example.TSP2;
import org.jenetics.internal.util.Equality;
import org.jenetics.internal.util.Hash;
import org.jenetics.util.MSeq;
import org.jenetics.util.RandomRegistry;

import java.util.Random;

import static java.lang.String.format;
import static org.jenetics.internal.math.random.indexes;

/**
 * The {@code SwapMutation} changes the order of genes in a chromosome, with the
 * hope of bringing related genes closer together, thereby facilitating the
 * production of building blocks. This mutation operator can also be used for
 * combinatorial problems, where no duplicated genes within a chromosome are
 * allowed, e.g. for the TSP.
 *
 * @author <a href="mailto:franz.wilhelmstoetter@gmx.at">Franz Wilhelmstötter</a>
 * @since 1.0
 * @version 3.0
 */
public class TSPMutator1C2<
	G extends Gene<?, EnumGene<Integer>>,
	C extends Comparable<? super C>
>
	extends Mutator<EnumGene<Integer>, C>
{

	/**
	 * Constructs an alterer with a given recombination probability.
	 *
	 * @param probability the crossover probability.
	 * @throws IllegalArgumentException if the {@code probability} is not in the
	 *          valid range of {@code [0, 1]}.
	 */
	public TSPMutator1C2(final double probability) {
		super(probability);
	}

	/**
	 * Default constructor, with default mutation probability
	 * ({@link AbstractAlterer#DEFAULT_ALTER_PROBABILITY}).
	 */
	public TSPMutator1C2() {
		this(DEFAULT_ALTER_PROBABILITY);
	}

	/*@Override
	public int alter(
			final Population<EnumGene<Integer>, C> population,
			final long generation
	) {return 0;}*/

/*	@Override
	public int alter(
			final Population<EnumGene<Integer>, C> population,
			final long generation
	) {

		assert(population != null) : "Not null is guaranteed from base class.";

		final double p = _probability;
		final IntRef alterations = new IntRef(0);

		indexes(RandomRegistry.getRandom(), population.size(), p).forEach(i -> {
			final Phenotype<EnumGene<Integer>, C> pt = population.get(i);

			final Genotype<EnumGene<Integer>> gt = pt.getGenotype();
			final Genotype<EnumGene<Integer>> mgt = mutate(gt, p, alterations);

			final Phenotype<EnumGene<Integer>, C> mpt = pt.newInstance(mgt, generation);
			population.set(i, mpt);
		});

		return alterations.value;
	}*/





	/**
	 * Swaps the genes in the given array, with the mutation probability of this
	 * mutation.
	 */
	//@Override
	/*protected int mutate2(final MSeq<EnumGene<Integer>> genes, final double p) {
		final Random random = RandomRegistry.getRandom();

		//Genotype<EnumGene<Integer>> g = new Genotype<EnumGene<Integer>>(genes.toISeq());

		//Chromosome<EnumGene<Integer>> c = Chromosome<EnumGene<Integer>>;
		//Genotype<EnumGene<Integer>> genotype = new Genotype<>(genes);


		return genes.length() > 1
			? (int)indexes(random, genes.length(), p)
				.peek(i -> genes.swap(i, random.nextInt(genes.length())))
				.count()
			: 0;
	}*/

	@Override
	protected int mutate(final MSeq<EnumGene<Integer>> genes, final double p) {
		final Random random = RandomRegistry.getRandom();

		//Chromosome<EnumGene<Integer>> c = Chromosome<EnumGene<Integer>>;
		//Genotype<EnumGene<Integer>> genotype = new Genotype<>(genes);

		System.out.println("ASD");
		System.out.println("test Dist2: " + TSP2.dist2(genes.toISeq()));





		return genes.length() > 1
				? (int)indexes(random, genes.length(), p)
				.peek(i -> genes.swap(i, random.nextInt(genes.length())))
				.count()
				: 0;
	}





	@Override
	public int hashCode() {
		return Hash.of(getClass()).and(super.hashCode()).value();
	}

	@Override
	public boolean equals(final Object obj) {
		return Equality.of(this, obj).test(super::equals);
	}

	@Override
	public String toString() {
		return format("%s[p=%f]", getClass().getSimpleName(), _probability);
	}

}
