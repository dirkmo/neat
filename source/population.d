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
        foreach( ref i; individuals ) { 
            i = new Individual( pool, true /*createConPhenotype*/ );
        }   
    }

    /// kill individuals with lowest fitness, fill up with
    /// offspring
    void selection() {
        //individuals.sort!( (a,b) => abs(a.fitness) < abs(b.fitness) )();
        auto species = new MedoidClassification!Individual(individuals, popsize / 10, float.nan);
        float dist = float.max;
        while( species.getMeanDistance() < dist ) {
            dist = species.getMeanDistance();
            species.doClustering();
        }

        uint biggestCluster, count;
        Individual[][] clusters;
        foreach( clidx; 0..species.getClusterCount() ) {
            uint len = cast(uint)species.getCluster(clidx).length;
            if( count < len ) {
                count = len;
                biggestCluster = clidx;
            }
            clusters ~= species.getCluster(clidx);
        }



        uint surplus = individuals.length - popsize;
        Individual[]
        while( surplus-- ) {

        }
        species.getCluster(clidx)

        Individual[][] newIndivduals;
        foreach( clidx; 0..species.getClusterCount() ) {
            auto ind = species.getCluster(clidx);
            ind.sort!( (a,b) => abs(a.fitness) < abs(b.fitness) )();
            uint survival = cast(uint)(ind.length * survival_rate);
            uint oldLength = cast(uint)ind.length;
            uint surplus;
            if( survival < 1 ) {
                newIndivduals[clidx] ~= ind[0];
                newIndivduals[clidx] ~= ind[0].crossOver( ind[0] );
                surplus++;
            } else {
                uint i;
                ind.length = survival;
                while(ind.length < oldLength) {
                    const uint p2 = uniform(0, survival);
                    auto offspring = ind[i++].crossOver(ind[p2]);
                    ind ~= offspring;
                }
                newIndivduals ~= ind;
            }
        }

        individuals = newIndivduals;

/*
        writefln("Best: %s, worst: %s, median: %s, average: %s", 
            individuals[0].fitness, individuals[$-1].fitness, individuals[$/2].fitness, average());
        uint survival = cast(uint)(individuals.length * survival_rate);
        uint oldLength = cast(uint)individuals.length;
        individuals.length = survival;

        uint i;
        while(individuals.length < oldLength ) {
            const uint p2 = uniform(0, survival);
            auto offspring = individuals[i++].crossOver(individuals[p2]);
            individuals ~= offspring;
        }
*/
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