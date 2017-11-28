module neat.population;

import neat.genepool;
import neat.individual;

import std.algorithm;
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

    void selection() {
        individuals.sort!( (a,b) => abs(a.fitness) < abs(b.fitness) )();
        uint survival = cast(uint)(individuals.length * survival_rate);
        uint newInd = cast(uint)individuals.length - survival;
        individuals.length = survival;
        foreach( i; 0..newInd ) {
            uint p2 = uniform(newInd, survival);
            auto offspring = individuals[i].crossOver(individuals[p2]);
            individuals ~= offspring;
        }
    }
        
//private:

    Genepool pool;
    Individual[] individuals;

    float survival_rate = 0.8;
}