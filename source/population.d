module neat.population;

import neat.genepool;
import neat.individual;

import std.algorithm;
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
        foreach( ref i; individuals ) { 
            i = new Individual( pool, true /*createConPhenotype*/ );
        }   
    }


    this( uint popsize, string topology, bool recurrent ) {
        pool = new Genepool(topology, recurrent);
        individuals.length = popsize;
        foreach( ref i; individuals ) { 
            i = new Individual( pool, true /*createConPhenotype*/ );
        }   
    }

    /// kill individuals with lowest fitness, fill up with
    /// offspring
    void selection() {
        individuals.sort!( (a,b) => abs(a.fitness) < abs(b.fitness) )();
        writefln("Best: %s, worst: %s, median: %s, average: %s", 
            individuals[0].fitness, individuals[$-1].fitness, individuals[$/2].fitness, average());
        uint survival = cast(uint)(individuals.length * survival_rate);
        uint oldLength = cast(uint)individuals.length;
        individuals.length = survival;
        uint i;
        while(individuals.length < oldLength ) {
            uint p2 = uniform(0, survival);
            auto offspring = individuals[i++].crossOver(individuals[p2]);
            individuals ~= offspring;
        }
    }

    void mutation() {
        foreach( i; individuals ) {
            i.mutateWeight( 0.1f, 1.0f );
            //i.mutateSplitUpConnection(0.01f);
            //i.mutateAddConnection(0.1f);
        }
    }

    float average() {
        float avg = 0;
        foreach(i; individuals) {
            avg += i.fitness;
        }
        return avg / individuals.length;
    }
        
//private:

    Genepool pool;
    Individual[] individuals;

    float survival_rate = 0.90;
}