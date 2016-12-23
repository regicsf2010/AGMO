package implementations;

import ga.Chromosome;
import ga.Constants;
import interfaces.Fitness;
import interfaces.SurvivorSelection;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;

import comparators.CrowdingDistanceComparator;
import comparators.FactoryComparator;
import comparators.ParetoFrontComparator;

public class NSGAII implements Fitness, SurvivorSelection {
	/**
	 * This method calculates pareto front from all chromosomes.
	 * This information is used to calculate all the chromosomes fitness.
	 * This method has overall complexity O(MN²), where M is the number of objectives
	 * 		and N is the population size.
	 * @param chromosomes All the chromosomes.
	 */
	static public void fastNonDominatedSort(Chromosome chromosomes[]) {
		Domination S[] = Domination.createDomination(chromosomes.length);
		ArrayList<Integer> F = new ArrayList<Integer>();
		
		for (int p = 0; p < chromosomes.length; p++) {
			for (int q = 0; q < chromosomes.length; q++) {
				if(p == q)
					continue;
				if(chromosomes[p].dominates(chromosomes[q]))
					S[p].dominated.add(q); // sol. p dominates sol. q
				else if(chromosomes[q].dominates(chromosomes[p]))
					S[p].dominationCount++; // sol. q dominates sol. p
			}
			if(S[p].dominationCount == 0) {
				chromosomes[p].setParetoFront(1);
				F.add(p);
			}
		}
		
		int i = 1; // i is the front counter
		while(!F.isEmpty()) {
			ArrayList<Integer> Q = new ArrayList<Integer>();
			for (Integer p : F) {
				for (Integer q : S[p].dominated) {
					S[q].dominationCount--;
					if(S[q].dominationCount == 0){
						chromosomes[q].setParetoFront(i + 1);
						Q.add(q);
					}
				}
			}			
			i++;
			F = Q;
		}		
	}

	private void calculateCrowdingDistances(Chromosome chromosomes[]) {
		for (int i = 0; i < chromosomes.length; i++)
			chromosomes[i].setCrowdingDistance(0);

		Collections.sort(Arrays.asList(chromosomes), new ParetoFrontComparator());
		int firstParetoFront = chromosomes[0].getParetoFront();
		int lastParetoFront = chromosomes[chromosomes.length - 1].getParetoFront();
		int firstPos = 0, size = 0;
		
		for (int i = firstParetoFront; i <= lastParetoFront; i++) {
			
			size = getParetoFrontSize(chromosomes, firstPos, i);
			Chromosome I[] = new Chromosome[size];
			System.arraycopy(chromosomes, firstPos, I, 0, I.length);
			
			if(size > 1) {
				for (int obj = 0; obj < Constants.NOBJECTIVES; obj++) {
					Collections.sort(Arrays.asList(I), FactoryComparator.getComparator(obj));
					
					I[0].setCrowdingDistance(Double.POSITIVE_INFINITY);
					I[I.length - 1].setCrowdingDistance(Double.POSITIVE_INFINITY);
					
					for (int j = 1; j < I.length - 1; j++) {
						I[j].setCrowdingDistance(I[j].getCrowdingDistance() + I[j + 1].getObjective(obj) - I[j - 1].getObjective(obj));
					}
				}
			
			} else {				
				I[0].setCrowdingDistance(Double.POSITIVE_INFINITY);
			}
			
			firstPos += size;			
		}		
	}	

	public static int getParetoFrontSize(Chromosome chromosomes[], int firstPos, int paretoFront) {
		int count = 0;
		for (int i = firstPos; i < chromosomes.length; i++) {
			if(chromosomes[i].getParetoFront() != paretoFront)
				break;
			count++;
		}		
		return count;
	}

	@Override
	public void calculateFitness(Chromosome chromosomes[]) {
		Chromosome.calculateAttributes(chromosomes);
		NSGAII.fastNonDominatedSort(chromosomes);
		this.calculateCrowdingDistances(chromosomes);
		// Finally calculates fitness of all chromosomes
		for (int i = 0; i < chromosomes.length; i++) {
			chromosomes[i].setFitness( (1f / chromosomes[i].getParetoFront()) + chromosomes[i].getCrowdingDistance() );
		}
	}

	@Override
	public void doSelection(Chromosome chromosomes[], Chromosome selectedChromosomes[]) {
		Chromosome all[] = new Chromosome[Constants.NPOPULATION * 2];
		for (int i = 0; i < Constants.NPOPULATION; i++) {
			all[i] = chromosomes[i];
			all[Constants.NPOPULATION + i] = selectedChromosomes[i];
		}
		
		NSGAII.fastNonDominatedSort(all);
		Collections.sort(Arrays.asList(all), new ParetoFrontComparator());
		
		Chromosome survivors[] = new Chromosome[Constants.NPOPULATION];
		
		int i = 1, firstPos = 0, size = 0, missingChromosomes = Constants.NPOPULATION;
		
		while(missingChromosomes != 0){
			size = getParetoFrontSize(all, firstPos, i);
			Chromosome F[] = new Chromosome[size];
			System.arraycopy(all, firstPos, F, 0, F.length);
			
			this.calculateCrowdingDistances(F);
			
			if(size <= missingChromosomes) {
				System.arraycopy(F, 0, survivors, firstPos, F.length);
				missingChromosomes -= size;
			} else {
				Collections.sort(Arrays.asList(F), new CrowdingDistanceComparator());
				System.arraycopy(F, 0, survivors, firstPos, missingChromosomes);
				missingChromosomes = 0;
			}			
			firstPos += size;
			i++;
		}
		
		System.arraycopy(survivors, 0, chromosomes, 0, Constants.NPOPULATION);
		// this is done because one solution may change its pareto front position or its crowding distance 
		this.calculateOnlyFitness(chromosomes);
	}
	
	private void calculateOnlyFitness(Chromosome chromosomes[]) {
		for (int i = 0; i < Constants.NPOPULATION; i++) 
			chromosomes[i].setFitness( (1f / chromosomes[i].getParetoFront()) + chromosomes[i].getCrowdingDistance() );	
	}

	private static class Domination {
		public int dominationCount;
		public ArrayList<Integer> dominated = new ArrayList<Integer>();
		
		public static Domination[] createDomination(int size){
			Domination d[] = new Domination[size];
			for (int i = 0; i < size; i++) 
				d[i] = new Domination();				
			return d;
		}
	}	
}
