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

import org.jenetics.example.RealFunction3;
import org.jenetics.internal.util.Equality;
import org.jenetics.internal.util.Hash;
import org.jenetics.internal.util.IntRef;
import org.jenetics.util.MSeq;
import org.jenetics.util.RandomRegistry;

import java.util.List;

import static java.lang.String.format;
import static org.jenetics.internal.math.random.indexes;

/**
 * This class is for mutating a chromosomes of an given population. There are
 * two distinct roles mutation plays
 * <ul>
 *     <li>Exploring the search space. By making small moves mutation allows a
 *     population to explore the search space. This exploration is often slow
 *     compared to crossover, but in problems where crossover is disruptive this
 *     can be an important way to explore the landscape.
 *     </li>
 *     <li>Maintaining diversity. Mutation prevents a population from
 *     correlating. Even if most of the search is being performed by crossover,
 *     mutation can be vital to provide the diversity which crossover needs.
 *     </li>
 * </ul>
 *
 * <p>
 * The mutation probability is the parameter that must be optimized. The optimal
 * value of the mutation rate depends on the role mutation plays. If mutation is
 * the only source of exploration (if there is no crossover) then the mutation
 * rate should be set so that a reasonable neighborhood of solutions is explored.
 * </p>
 * The mutation probability <i>P(m)</i> is the probability that a specific gene
 * over the whole population is mutated. The number of available genes of an
 * population is
 * <p>
 * <img src="doc-files/mutator-N_G.gif" alt="N_P N_{g}=N_P \sum_{i=0}^{N_{G}-1}N_{C[i]}" >
 * </p>
 * where <i>N<sub>P</sub></i>  is the population size, <i>N<sub>g</sub></i> the
 * number of genes of a genotype. So the (average) number of genes
 * mutated by the mutation is
 * <p>
 * <img src="doc-files/mutator-mean_m.gif" alt="\hat{\mu}=N_{P}N_{g}\cdot P(m)" >
 * </p>
 *
 * @author <a href="mailto:franz.wilhelmstoetter@gmx.at">Franz Wilhelmstötter</a>
 * @since 1.0
 * @version 3.0
 */
public class Mutator_M6<
	G extends Gene<Double, DoubleGene>,
	C extends Comparable<? super C>
>
	extends AbstractAlterer<DoubleGene, C>
{

	/**
	 * Construct a Mutation object which a given mutation probability.
	 *
	 * @param probability Mutation probability. The given probability is
	 *         divided by the number of chromosomes of the genotype to form
	 *         the concrete mutation probability.
	 * @throws IllegalArgumentException if the {@code probability} is not in the
	 *          valid range of {@code [0, 1]}..
	 */
	public Mutator_M6(final double probability) {
		super(probability);
	}

	/**
	 * Default constructor, with probability = 0.01.
	 */
	public Mutator_M6() {
		this(0.01);
	}

	/**
	 * Concrete implementation of the alter method.
	 */
	@Override
	public int alter(
		final Population<DoubleGene, C> population,
		final long generation
	) {
		assert(population != null) : "Not null is guaranteed from base class.";

		//final double p = pow(_probability, 1.0/3.0);
		final double p = _probability;

		final IntRef alterations = new IntRef(0);
		//System.out.println("prob:"+_probability);
		//System.out.println("p:"+p);

		indexes(RandomRegistry.getRandom(), population.size(), p).forEach(i -> {
			final Phenotype<DoubleGene, C> pt = population.get(i);

			final Genotype<DoubleGene> gt = pt.getGenotype();
			final Genotype<DoubleGene> mgt = mutate(gt, p, alterations);

			final Phenotype<DoubleGene, C> mpt = pt.newInstance(mgt, generation);
			population.set(i, mpt);
		});


		return alterations.value;
	}

	private Genotype<DoubleGene> mutate(
		final Genotype<DoubleGene> genotype,
		final double p,
		final IntRef alterations
	) {
		final MSeq<Chromosome<DoubleGene>> chromosomes = genotype.toSeq().copy();


		alterations.value +=
			indexes(RandomRegistry.getRandom(), genotype.length(), p)
				.map(i -> mutate(chromosomes, i, p))
				.sum();

		return genotype.newInstance(chromosomes.toISeq());
	}

	private int mutate(final MSeq<Chromosome<DoubleGene>> c, final int i, final double p) {

		double change = 0.01;

		System.out.println("i: "+i);
		final Chromosome<DoubleGene> chromosome = c.get(i);
		MSeq<DoubleGene> genes = chromosome.toSeq().copy();

		final MSeq<Chromosome<DoubleGene>> cM = c.copy();
		final MSeq<Chromosome<DoubleGene>> cW = c.copy();

		//c.copy();

		final MSeq<DoubleGene> genesM = chromosome.toSeq().copy();
		final MSeq<DoubleGene> genesW = chromosome.toSeq().copy();

		double val = genes.get(0).doubleValue();

		/*System.out.println("val:"+val);

		System.out.println("genesM: "+genesM.get(0).doubleValue());
		System.out.println("genesW: "+genesW.get(0).doubleValue());

		System.out.println("Test0:");
		System.out.println(evalG(c));
		System.out.println(evalG(cM));
		System.out.println(evalG(cW));*/

		genesM.set(0, DoubleGene.of(Math.max(val-change, RealFunction3.getMinV()), RealFunction3.getMinV(), RealFunction3.getMaxV()));
		genesW.set(0, DoubleGene.of(Math.min(val+change, RealFunction3.getMaxV()), RealFunction3.getMinV(), RealFunction3.getMaxV()));

		//System.out.println("genesM: " + genesM.get(0).doubleValue());
		//System.out.println("genesW: " + genesW.get(0).doubleValue());
		
		//cM.set(i, cM.get(i).newInstance(genesM.toISeq()));
		//cW.set(i, cW.get(i).newInstance(genesW.toISeq()));

		//cM.set(i, new DoubleChromosome(genesM.toISeq()));
		//cW.set(i, new DoubleChromosome(genesW.toISeq()));

		cM.set(i, DoubleChromosome.of(genesM.get(0)));
		cW.set(i, DoubleChromosome.of(genesW.get(0)));


		System.out.println("Test1:");
		System.out.println(evalG(c));
		System.out.println(evalG(cM));
		System.out.println(evalG(cW));

		/*
		if(genesM==genesW){
			System.out.println("Test2 F");
		}else{
			System.out.println("Test2 OK");
		}

		cM.set(0, DoubleChromosome.of(DoubleGene.of(5, -5, 5)));
		cM.set(1, DoubleChromosome.of(DoubleGene.of(5, -5, 5)));

		System.out.println("Test3:");
		System.out.println(evalG(c));
		System.out.println(evalG(cM));
		System.out.println(evalG(cW));*/


		int mutations = 0;

//		/*final int mutations = mutate(genes, p);
//		if (mutations > 0) {
//			c.set(i, chromosome.newInstance(genes.toISeq()));
//		}
//
//		return mutations;*/

		int option = 0;
		double currentMin = evalG(c);
		if(evalG(cM)<currentMin){
			option = 1;
			mutations++;
		}
		if(evalG(cW)<currentMin){
			option = 2;
			mutations++;
		}

		switch(option){
			case 0: c.set(i, chromosome.newInstance(genes.toISeq())); break;
			case 1: c.set(i, chromosome.newInstance(genesM.toISeq())); break;
			case 2: c.set(i, chromosome.newInstance(genesW.toISeq())); break;
		}

		System.out.println("mutations: "+mutations);

		return mutations;

	}

	/**
	 * <p>
	 * Template method which gives an (re)implementation of the mutation class
	 * the possibility to perform its own mutation operation, based on a
	 * writable gene array and the gene mutation probability <i>p</i>.
	 *
	 * @param genes the genes to mutate.
	 * @param p the gene mutation probability.
	 * @return the number of performed mutations
	 */
	protected int mutate(final MSeq<DoubleGene> genes, final double p) {
		//Double zxc = new Double(5);

		//System.out.println("p: " +p);

		//final Number zxc = 1.5;
		//final Double zxc2 = 1.5;
		//final Number zxc3 = new Double(2.5);
		//final Double zxc4 = -5.0;



		int ret = (int)indexes(RandomRegistry.getRandom(), genes.length(), p)
			.peek(i -> genes.set(i, genes.get(i).newInstance()))
			//.peek(i -> genes.set(i, genes.get(i).newInstance(zxc4)))
				//.peek(i -> fooHelper2(genes))
				//.peek(i -> foo3(genes))
			.count();
		//System.out.println("ret:"+ret);
		return ret;



//		return (int)indexes(RandomRegistry.getRandom(), genes.length(), p)
//			//.peek(i -> genes.set(i, genes.get(i).newInstance()))
//			.peek(i -> genes.set(i, genes.get(i).newInstance(zxc4)))
//				//.peek(i -> fooHelper2(genes))
//				//.peek(i -> foo3(genes))
//				.count();

		//genes.set(0, genes.get(0).newInstance(zxc4));
		//genes.set(1, genes.get(1).newInstance(zxc4));
		//return 1;

	    //System.out.println(x.doubleValue() + ", " + y.doubleValue());
		//genes.set(1, x);

		//DoubleGene x = genes.get(0);
		//DoubleGene y = genes.get(1);

		//genes.set(1, upgrade1(x,y));

		//return 1;
	}

	/*private DoubleGene upgrade1(DoubleGene x, DoubleGene uy){
		double dx = x.doubleValue();
		double v = uy.doubleValue();

		//System.out.println("v:"+v);
		double vW = testD(v + 0.005);
		double vM = testD(v - 0.005);


		if (RealFunction3.eval2(dx, vW) < RealFunction3.eval2(dx, v)){
		v = vW;
		}
		if (RealFunction3.eval2(dx, vM) < RealFunction3.eval2(dx, v)){
			v = vM;
		}
		//System.out.println("v2:"+v);

		return DoubleGene.of(v,RealFunction3.getMinV(),RealFunction3.getMaxV());
	}*/

	private double testD(double a){
	if(a < RealFunction3.getMinV()) return RealFunction3.getMinV();
	if(a > RealFunction3.getMaxV()) return RealFunction3.getMaxV();
	return a;
	}

	private double evalG(MSeq<Chromosome<DoubleGene>> c){
		Genotype<DoubleGene> g = new Genotype<DoubleGene>(c.toISeq());
		//g.newInstance(c.toISeq());
		return RealFunction3.eval(g);
	}


	void foo3(MSeq<?> i) {
		fooHelper2(i);
	}


	void foo(List<?> i) {
		fooHelper(i);
	}

	private <T> void fooHelper(List<T> asd){
		asd.set(0, asd.get(0));
	}
/*
	private <T extends DoubleGene> void fooHelper2(int ind, final MSeq<T> gen){
		Number asd = 1.5;
		gen.set(ind, gen.get(ind).newInstance());
	}
*/
    void foo2(MSeq<? extends DoubleGene> i) {
	fooHelper2(i);
}

	private <G> void fooHelper2(final MSeq<G> asd){
		//System.out.println(asd.get(0).getClass());
		asd.set(0, asd.get(0));

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
