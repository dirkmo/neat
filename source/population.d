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
        classificator = new SpeciesClassificator(cast(Phenotype[])individuals, 2.0f);
    }


    this( uint popsize, string topology, bool recurrent ) {
        pool = new Genepool(topology, recurrent);
        individuals.length = popsize;
        this.popsize = popsize;        
        foreach( ref i; individuals ) { 
            i = new Individual( pool, true /*createConPhenotype*/ );
        }
        classificator = new SpeciesClassificator(cast(Phenotype[])individuals, 2.0f);
    }

    /// kill individuals with lowest fitness, fill up with offspring
    void selection() {
        Individual[] newIndividuals;
        float totalFitness = classificator.totalFitness(individuals);
        writefln("Totalfitness: %s", totalFitness);
        float[] speciesFitness;

        // remove weakest members from species and fill up with new offspring.
        // Member count of a species is dependend on shared fitness.
        foreach( speciesIdx; 0..classificator.numberOfSpecies ) {
            // get all members of this species
            auto speciesMembers = individuals.filter!(i => i.species == speciesIdx).array;
            // sort by fitness
            speciesMembers.sort!( (a, b) => a.fitness < b.fitness );
            // get fitness of species
            speciesFitness ~= classificator.speciesFitness(individuals, speciesIdx);
            writefln("SpeciesFitness %s: %s", speciesIdx, speciesFitness[$-1]);
            // save members count of this species
            const uint oldlen = cast(uint)speciesMembers.length;
            // change species member count according to species fitness in relation to population fitness
            speciesMembers.length = cast(uint)(popsize * speciesFitness[speciesIdx] / totalFitness);
            // create offspring from parent1 and parent2
            uint parent1, parent2;
            uint survival;
            if( oldlen >= speciesMembers.length ) {
                // member count does not grow
                survival = cast(uint)(speciesMembers.length * survival_rate);
            } else {
                // member count grows
                survival = cast(uint)(oldlen * survival_rate);
            }
            // [0..survival] survive, [survival..speciesMembers.length) replaced by new offspring
            for( uint i = survival; i < speciesMembers.length; i++ ) {
                if( survival == 0 ) {
                    hier weiter
                    
                } else {
                    // parent2 is randomly picked
                    parent2 = uniform(0, survival);
                }
                speciesMembers[i] = speciesMembers[parent1].crossOver(speciesMembers[parent2]);
                classificator.determineSpecies(speciesMembers[i], individuals);
                parent1++; // parent1 just counting up
            }
            newIndividuals ~= speciesMembers;
        }
        individuals = newIndividuals;
        writeln("Individual count: ", individuals.length);
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