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

        foreach( sp; speciesClassificator.range() ) {
            const uint spIdx = sp.index;
            // get members of species
            auto members = individuals.filter!(a => a.species == spIdx).array;
            // fittest individuals first
            members.sort!( (a,b) => a.fitness > b.fitness );
            // kill the unfit
            uint survival = cast(uint)(members.length * survivalRate);
            if( survival < 1 ) {
                // kill species
                speciesClassificator.extinctSpecies(spIdx);
                writefln("Species %s goes extinct.", spIdx);
            } else {
                members.length = survival;
                uint parent1;
                while( members.length < sp.nextGenMemberCount ) {
                    uint parent2 = uniform(0, survival);
                    auto offspring = members[parent1].crossOver(members[parent2]);
                    members ~= offspring;
                    parent1 = (parent1 + 1) % survival;
                }
                newIndividuals ~= members;
            }
        }

        individuals = newIndividuals;
        speciesClassificator.reassign(individuals);
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

    enum survivalRate = 0.90f;
    enum SpeciesThreshold = 2.0f;
}