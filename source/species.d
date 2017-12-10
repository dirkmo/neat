module neat.species;

import neat.phenotype;

import std.algorithm;
import std.random;
import std.range;
import std.stdio;
import std.typecons;

class SpeciesClassificator {
    this( Phenotype[]* individuals, float thresh) {
        this.individuals = individuals;
        this.thresh = thresh;
        prototypes ~= (*individuals)[0];
        foreach(i; *individuals) {
            auto best = bestMatch(i);
            if( best[0] > thresh ) {
                prototypes ~= i.clone();
                i.species = cast(uint)prototypes.length - 1;
                writeln("Added species ", i.species);
            } else {
                i.species = best[1];
            }
        }
    }


    private auto bestMatch( Phenotype ind ) {
        float bestDist = float.max;
        uint bestIdx = uint.max;
        foreach(idx, p; prototypes) {
            float dist = ind.distance(ind, 1.0f, 1.0f, 0.5f );
            if( dist < bestDist ) {
                bestDist = dist;
                bestIdx = cast(uint)idx;
            }
        }
        return tuple(bestDist, bestIdx);
    }

    void updatePrototypes(Phenotype[] individuals) {
        foreach( idx, ref proto; prototypes ) {
            auto r = individuals.find!( i => i.species == idx )().array;
            proto = r[uniform(0,$)];
        }
    }


    float sharedFitness(uint species) {
        auto members = (*individuals).filter!(i=>i.species == species)();
        uint count;
        float fitness;
        foreach(m;members) {
            fitness += m.fitness;
            count++;
        }
        return fitness / count;
    }

    Phenotype[]* individuals;
    Phenotype[] prototypes;
    float thresh;
}