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
import std.range;

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
        speciesClassificator = new SpeciesClassificator(individuals, SpeciesThreshold, 5);
    }


    this( uint popsize, string topology, bool recurrent ) {
        pool = new Genepool(topology, recurrent);
        individuals.length = popsize;
        this.popsize = popsize;        
        foreach( ref i; individuals ) { 
            i = new Individual( pool, true /*createConPhenotype*/ );
        }
        speciesClassificator = new SpeciesClassificator(individuals, SpeciesThreshold, 5);
    }

    /// kill individuals with lowest fitness and fill up empty spaces with
    /// new offspring
    void selection() {
        // make sure that all individuals belong to a species
        speciesClassificator.update(individuals);
        // calculate fitness of species
        speciesClassificator.calculateFitness();
        // calculate size of species for the next generation
        speciesClassificator.calculateNextGenSpeciesSize(popsize);

        Individual[] newIndividuals;
        
        //writeln("Starting selection");
        foreach( sp; speciesClassificator.range() ) {
            const uint spIdx = sp.index;
            //writefln("\n***Species %s***", spIdx);
            // get members of species
            auto members = individuals.filter!(a => a.species == spIdx).array;

            assert(members.length == sp.memberCount);

            // fittest individuals first
            members.sort!( (a,b) => a.fitness > b.fitness );

            uint survival = cast(uint)(sp.nextGenMemberCount * survivalRate);
            if( sp.isNew ) {
                // new / immature species, won't extinct
                if( survival == 0 ) survival = 1;
            }
            if( sp.memberCount < survival ) survival = sp.memberCount;
            if( survival > 0 ) {
                // kill the unworthy
                members.length = survival;
                // create offspring
                uint parent1;
                while( members.length < sp.nextGenMemberCount ) {
                    uint parent2 = uniform(0, survival);
                    if( members[parent2] is null ) {
                        writefln("len: %s, survival: %s, parent2: %s, sp.memberCount: %s",
                                 members.length, survival, parent2, sp.memberCount);
                        throw new Exception("warum nur???");
                    }
                    auto offspring = members[parent1].crossOver(members[parent2]);
                    members ~= offspring;
                    parent1 = (parent1 + 1) % survival;
                }
                newIndividuals ~= members;
            }
        }

        individuals = newIndividuals;
        // pick new prototypes
        speciesClassificator.pickNewPrototypes(individuals);
        speciesClassificator.age();
        //writeln("Individual count: ", individuals.length);
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

    Individual best() {
        float max = 0;
        ulong idx;
        foreach(lidx, i; individuals) {
            if( i.fitness > max ) {
                max = i.fitness;
                idx = lidx;
            }
        }
        return individuals[idx];
    }


//private:

    Genepool pool;
    Individual[] individuals;
    uint popsize;

    SpeciesClassificator speciesClassificator;

    enum survivalRate = 0.90f;
    enum SpeciesThreshold = 3.0f;
}