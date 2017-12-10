module neat.population;

import neat.genepool;
import neat.individual;
import neat.phenotype;
import neat.species;

import std.algorithm;
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
        classificator = new SpeciesClassificator();
    }


    this( uint popsize, string topology, bool recurrent ) {
        pool = new Genepool(topology, recurrent);
        individuals.length = popsize;
        this.popsize = popsize;        
        foreach( ref i; individuals ) { 
            i = new Individual( pool, true /*createConPhenotype*/ );
        }
        classificator = new SpeciesClassificator();
    }

    /// kill individuals with lowest fitness, fill up with offspring
    void selection() {
        individuals.sort!( (a,b) => abs(a.fitness) < abs(b.fitness) )();
        uint species = classificator.assignSpecies(cast(Phenotype[])individuals, 10.0f);
        foreach( s; 0..species ) {
            auto members = individuals.filter!(i=>i.species == s)();
            
        }
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