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
            true /*fullyConnected*/,
            recurrent /*enableRecurrentNets*/
        );  
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
        uint newInd = cast(uint)individuals.length - survival;
        individuals.length = survival;
        foreach( i; 0..newInd ) {
            uint p2 = uniform(newInd, survival);
            auto offspring = individuals[i].crossOver(individuals[p2]);
            individuals ~= offspring;
        }
    }

    void mutation() {
        foreach( i; individuals ) {
            i.mutateWeight( 0.2f, 1.5f );
            i.mutateSplitUpConnection(0.1f);
            i.mutateAddConnection(0.1f);
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

    float survival_rate = 0.9;
}