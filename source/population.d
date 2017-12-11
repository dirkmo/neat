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
        classificator = new SpeciesClassificator(cast(Phenotype[])individuals, 10.0f);
    }


    this( uint popsize, string topology, bool recurrent ) {
        pool = new Genepool(topology, recurrent);
        individuals.length = popsize;
        this.popsize = popsize;        
        foreach( ref i; individuals ) { 
            i = new Individual( pool, true /*createConPhenotype*/ );
        }
        classificator = new SpeciesClassificator(cast(Phenotype[])individuals, 10.0f);
    }

    /// kill individuals with lowest fitness, fill up with offspring
    void selection() {
        Individual[] newIndividuals;
        float totalFitness = classificator.totalFitness(individuals);
        float[] speciesFitness;

        // remove weakest members from species and fill up with new offspring.
        // Member count of a species is dependend of shared fitness.
        foreach( speciesIdx; 0..classificator.numberOfSpecies ) {
            // get all members of this species
            auto speciesMembers = individuals.filter!(i => i.species == speciesIdx).array;
            speciesMembers.sort!( (a, b) => a.fitness < b.fitness );
            speciesFitness[speciesIdx] = classificator.speciesFitness(individuals, speciesIdx);
            const uint oldlen = cast(uint)speciesMembers.length;
            // change species member count according to fitnesss
            speciesMembers.length = cast(uint)(popsize * speciesFitness[speciesIdx] / totalFitness);

            if( speciesMembers.length < oldlen ) {
                
            } else {

            }

            newIndividuals ~= speciesMembers;
        }

        individuals = newIndividuals;
        classificator.updatePrototypes(cast(Phenotype[])individuals);
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
    SpeciesClassificator classificator;
    uint popsize;

    float survival_rate = 0.90;
}