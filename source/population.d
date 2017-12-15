module neat.population;

import neat.genepool;
import neat.individual;
import neat.phenotype;
import neat.species;

import std.algorithm;
import std.array;
import std.container;
import std.stdio;
import std.math;
import std.random;

class Population {                                                           
    this( uint popsize, uint inputs, uint outputs, bool recurrent ) { 
        pool = new Genepool(
            inputs, outputs,
            true /*addBias*/,
            true /*fullyConnected*/,
            recurrent /*enableRecurrentNets*/
        );
        individuals.length = popsize;
        this.popsize = popsize;
        foreach( ref i; individuals ) { 
            i = new Individual( pool, true /*createConPhenotype*/ );
        }
        speciesClassificator = new SpeciesClassificator(individuals, SpeciesThreshold);
    }


    this( uint popsize, string topology, bool recurrent ) {
        pool = new Genepool(topology, recurrent);
        individuals.length = popsize;
        this.popsize = popsize;        
        foreach( ref i; individuals ) { 
            i = new Individual( pool, true /*createConPhenotype*/ );
        }
        speciesClassificator = new SpeciesClassificator(individuals, SpeciesThreshold);
    }

    /// kill individuals with lowest fitness and fill up empty spaces with
    /// new offspring
    void selection() {
        // check that all individuals belong to a species
        speciesClassificator.update(individuals);
        // calculate fitness of species
        speciesClassificator.calculateFitness();
        // calculate size of species for the next generation
        speciesClassificator.calculateNextGenSpeciesSize(popsize);

        Individual[] newIndividuals;

        individuals = newIndividuals;
        speciesClassificator.reassign();
        writeln("Individual count: ", individuals.length);
    }

    void mutation() {
        foreach( i; individuals ) {
            i.mutateWeight( 0.1f, 1.0f );
            i.mutateSplitUpConnection(0.001f);
            i.mutateAddConnection(0.001f);
        }
    }

    float average() {
        float avg = 0;
        foreach(i; individuals) {
            avg += i.fitness;
        }
        return avg / individuals.length;
    }

    Individual first() {
        return individuals[0];
    }


//private:

    Genepool pool;
    Individual[] individuals;
    uint popsize;

    SpeciesClassificator speciesClassificator;

    enum survival_rate = 0.90f;
    enum SpeciesThreshold = 2.0f;
}