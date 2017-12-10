module neat.population;

import neat.genepool;
import neat.individual;
import neat.medoids;

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
    }


    this( uint popsize, string topology, bool recurrent ) {
        pool = new Genepool(topology, recurrent);
        individuals.length = popsize;
        this.popsize = popsize;        
        foreach( ref i; individuals ) { 
            i = new Individual( pool, true /*createConPhenotype*/ );
        }   
    }

    /// kill individuals with lowest fitness, fill up with
    /// offspring
    void selection() {
        //individuals.sort!( (a,b) => abs(a.fitness) < abs(b.fitness) )();

        auto clusterAlgorithm = new MedoidClassification!Individual(individuals, popsize / 10, float.nan);
        auto species = clusterAlgorithm.clusterAll();
        //uint biggestSpeciesIdx = cast(uint)species.maxIndex!( (a,b) => a.length < b.length )();
        /*

        uint surplus;

        Individual[][] newIndividuals;
        foreach(idx, sp; species) {
            Individual[] newInd;
            sp.sort!( (a,b) => abs(a.fitness) < abs(b.fitness) )();
            uint survival = cast(uint)(sp.length * survival_rate);
            uint oldLength = cast(uint)sp.length;
            if( survival < 1 ) {
                newInd ~= sp[0];
                newInd ~= sp[0].crossOver( sp[0] ); // clone
                surplus++;
            } else {
                // copy fittest
                newInd.copy(sp[0..survival]);                
                uint i;
                while(newInd.length < oldLength) {
                    const uint p2 = uniform(0, survival);
                    auto offspring = sp[i++].crossOver( sp[p2] );
                    newInd ~= offspring;
                }
            }
            newIndividuals ~= newInd;
        }
        individuals.length = 0;
        foreach(ni; newIndividuals) {
            individuals ~= ni;
        }
        writeln("Count: ", individuals.length);
*/

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
    
    void cluster(float dist) {

    }

//private:

    Genepool pool;
    Individual[] individuals;
    uint popsize;

    float survival_rate = 0.90;
}