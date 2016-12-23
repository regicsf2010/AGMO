package implementations;

import ga.Chromosome;
import ga.Constants;
import interfaces.Fitness;
import interfaces.SurvivorSelection;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;

import comparators.FitnessAscendantComparator;
import comparators.SpeaIIDistanceComparator;

public class SPEAII implements Fitness, SurvivorSelection {

	// This archive is static for convenient.
	public static ArrayList<Chromosome> archive = new ArrayList<>(); 
	
	@Override
	public void doSelection(Chromosome[] chromosomes, Chromosome[] selectedChromosomes) {
		// if it is in the first generation, archive is empty
		if(archive.isEmpty())
			archive = this.getAllNonDominated(chromosomes, selectedChromosomes);
		else { 
			// other generation, archive is not empty and 
			// the members have their fitness updated. 
			archive = this.getAllNonDominated(archive.toArray(new Chromosome[archive.size()]), selectedChromosomes);
		}
			
		if(archive.size() > Constants.NARCHIVE) {
			this.trucationArchive();
		} else {
			this.completeArchive(chromosomes, selectedChromosomes);			
		}

		for (int i = 0; i < archive.size(); i++) 
			chromosomes[i] = archive.get(i);
		
	}

	private void trucationArchive() {
		int k = (int) Math.sqrt(archive.size());
		
		while (archive.size() > Constants.NARCHIVE) {
			SpeaIIChromosomeDistance allDistances[] = SpeaIIChromosomeDistance.createSpeaIIChromosomeDistace(archive.size());
			
			for (int p = 0; p < archive.size(); p++) {
				double distances[] = new double[archive.size() - 1]; // p doesn't have distance to its self
				int count = 0;
				for (int q = 0; q < archive.size(); q++) {
					if(p != q) {
						distances[count++] = this.calculateDistance(archive.get(p), archive.get(q));
					}
				}	
				Arrays.sort(distances);
				allDistances[p].distance = distances[k];
				allDistances[p].chromosomePosition = p;
			}
			
			Collections.sort(Arrays.asList(allDistances), new SpeaIIDistanceComparator());

			/**
			 * Assumes position one has the smallest distance.
			 * But position two, three and so on can have the same distance of the first one. In this case,
			 * removes the second smallest distance. 
			 */
			int position = allDistances[0].chromosomePosition;
			double shortestDistance = allDistances[0].distance;			
			if(shortestDistance != allDistances[1].distance) {
				archive.remove(position);
				continue;
			}
			/**
			 * In this case, first and second position have the same distance.
			 * So, removes the second smallest distance (or the third smallest distance and so on).
			 */
			for (int i = 2; i < allDistances.length; i++) {
				if(shortestDistance != allDistances[i].distance) {
					position = allDistances[i].chromosomePosition;
					shortestDistance = allDistances[i].distance;	
					if(i + 1 < allDistances.length)
						if(shortestDistance != allDistances[i + 1].distance)
							break;						
				}				
			}			
			archive.remove(position);
		}
	}
	
	private void completeArchive(Chromosome[] chromosomes, Chromosome[] selectedChromosomes) {
		ArrayList<Chromosome> all = new ArrayList<>();

		for (int i = 0; i < Constants.NPOPULATION; i++) {
			all.add(chromosomes[i]);
			all.add(selectedChromosomes[i]);
		}
		
		Collections.sort(all, new FitnessAscendantComparator());
		int n = Constants.NARCHIVE - archive.size();
		
		
		for (int i = 0; i < all.size(); i++) {
			if(all.get(i).getRawFitness() != 0) {
				archive.add(all.get(i));
				n--;
			}
			if(n == 0)
				break;
		}
		
	}
	
	@Override
	public void calculateFitness(Chromosome[] chromosomes) {
		Chromosome.calculateAttributes(chromosomes);
		this.calculateStrength(chromosomes);
		this.calculateRawFitness(chromosomes);		
		this.calculateDensity(chromosomes);
		
		for (int i = 0; i < Constants.NPOPULATION; i++) 
			chromosomes[i].setFitness(chromosomes[i].getRawFitness() + chromosomes[i].getDensity());
		if(!archive.isEmpty())
			for (int i = 0; i < Constants.NARCHIVE; i++) {
				archive.get(i).setFitness(archive.get(i).getRawFitness() + archive.get(i).getDensity());
			}		
	}
	
	/**
	 * This method calculate the strength of population and archive individuals.
	 * The calculation is stored in the structures SpeaIIChromosomeInfo. 
	 * @param chromosomes The population individuals.
	 * @param speaIIChromosome Structure for storing the chromosomes strength.
	 * @param speaIIArchive Structure for storing the archive strength.
	 */
	private void calculateStrength(Chromosome[] chromosomes) {
		
		/**
		 * Calculates the chromosomes strengths.
		 */
		for (int p = 0; p < Constants.NPOPULATION; p++) {
			for (int q = 0; q < Constants.NPOPULATION; q++) {
				if(p != q) {
					if (chromosomes[p].dominates(chromosomes[q]))
						chromosomes[p].setStrength(chromosomes[p].getStrength() + 1);
				}				
			}
			if(!archive.isEmpty())
				
				for (int w = 0; w < Constants.NARCHIVE; w++) {
					if (chromosomes[p].dominates(archive.get(w)))
						chromosomes[p].setStrength(chromosomes[p].getStrength() + 1);
				}
		}
		
		/**
		 * Calculates the archive strengths.
		 */
		if(!archive.isEmpty())
			for (int p = 0; p < Constants.NARCHIVE; p++) {
				for (int q = 0; q < Constants.NARCHIVE; q++) {
					if (p != q) {
						if (archive.get(p).dominates(archive.get(q)))
							archive.get(p).setStrength(archive.get(p).getStrength() + 1);
					}
				}
				for (int w = 0; w < Constants.NPOPULATION; w++) {
					if (archive.get(p).dominates(chromosomes[w]))
						archive.get(p).setStrength(archive.get(p).getStrength() + 1);
				}
			}		
	}
	
	/**
	 * This method calculate the raw fitness of population and archive individuals.
	 * The calculation is stored in the structure SpeaIIChromosomeInfo. 
	 * @param chromosomes
	 * @param speaII
	 */
	private void calculateRawFitness(Chromosome[] chromosomes) {
		/**
		 * Calculates the raw fitness of all chromosomes.
		 */
		for (int p = 0; p < Constants.NPOPULATION; p++) {
			for (int q = 0; q < Constants.NPOPULATION; q++) {
				if(p != q) {
					if(chromosomes[q].dominates(chromosomes[p]))
						chromosomes[p].setRawFitness(chromosomes[p].getRawFitness() + chromosomes[q].getStrength());						
				}			
			}	
			if(!archive.isEmpty())
				for (int w = 0; w < Constants.NARCHIVE; w++) {
					if (archive.get(w).dominates(chromosomes[p]))
						chromosomes[p].setRawFitness(chromosomes[p].getRawFitness() + archive.get(w).getStrength());						
				}
		}
		
		/**
		 * Calculates the raw fitness of all archive members.
		 */		
		if(!archive.isEmpty())
			for (int p = 0; p < Constants.NARCHIVE; p++) {
				for (int q = 0; q < Constants.NARCHIVE; q++) {
					if (p != q) {
						if (archive.get(q).dominates(archive.get(p)))
							archive.get(p).setRawFitness(archive.get(p).getRawFitness() + archive.get(q).getStrength());							
					}
				}
				for (int w = 0; w < Constants.NPOPULATION; w++) {
					if (chromosomes[w].dominates(archive.get(p)))
						archive.get(p).setRawFitness(archive.get(p).getRawFitness() + chromosomes[w].getStrength());						
				}
			}
		
	}
	
	private void calculateDensity(Chromosome[] chromosomes) {
		int k = archive.isEmpty()? (int) Math.sqrt(Constants.NPOPULATION) : (int) Math.sqrt(Constants.NPOPULATION + Constants.NARCHIVE);
		
		for (int p = 0; p < Constants.NPOPULATION; p++) {
			double distances[] = archive.isEmpty()? new double[Constants.NPOPULATION - 1] : new double[Constants.NPOPULATION + Constants.NARCHIVE - 1]; // p doesn't have distance to its self 
			int count = 0;
	
			for (int q = 0; q < Constants.NPOPULATION; q++) {
				if(p != q) {
					distances[count++] = this.calculateDistance(chromosomes[p], chromosomes[q]);
				}				
			}
			
			if(!archive.isEmpty())
				for (int w = 0; w < Constants.NARCHIVE; w++) {
					distances[count++] = this.calculateDistance(chromosomes[p], archive.get(w));
				}
			
			Arrays.sort(distances);
			chromosomes[p].setDensity(1 / (2 + distances[k]));
		}
		
		if(!archive.isEmpty())
			for (int p = 0; p < Constants.NARCHIVE; p++) {
				double distances[] = new double[Constants.NPOPULATION + Constants.NARCHIVE - 1]; // p doesn't have distance to its self
				int count = 0;
				
				for (int q = 0; q < Constants.NARCHIVE; q++) {
					if (p != q) {
						distances[count++] = this.calculateDistance(archive.get(p), archive.get(q));
					}
				}
				
				for (int w = 0; w < Constants.NPOPULATION; w++) {
					distances[count++] = this.calculateDistance(archive.get(p), chromosomes[w]);
				}
				
				Arrays.sort(distances);
				archive.get(p).setDensity(1 / (2 + distances[k]));
			}
	}
	
	private double calculateDistance(Chromosome c1, Chromosome c2) {
		return Math.sqrt(Math.pow(c1.getObjective(0) - c2.getObjective(0), 2) + Math.pow(c1.getObjective(1) - c2.getObjective(1), 2));
	}
	
	private ArrayList<Chromosome> getAllNonDominated(Chromosome chromosomes[], Chromosome selectedChromosomes[]) {
		ArrayList<Chromosome> archiveL = new ArrayList<>();
		for (int p = 0; p < Constants.NPOPULATION; p++) {
			if(chromosomes[p].getRawFitness() == 0)
				archiveL.add(chromosomes[p]);
			if(selectedChromosomes[p].getRawFitness() == 0)
				archiveL.add(chromosomes[p]);			
		}				
		return archiveL;		
	}
	
	public static class SpeaIIChromosomeDistance{
		public int chromosomePosition;
		public double distance;		
		
		public static SpeaIIChromosomeDistance[] createSpeaIIChromosomeDistace(int size){
			SpeaIIChromosomeDistance d[] = new SpeaIIChromosomeDistance[size];
			for (int i = 0; i < size; i++) 
				d[i] = new SpeaIIChromosomeDistance();				
			return d;
		}
	}
}
