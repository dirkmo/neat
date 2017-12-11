module neat.species;

import neat.phenotype;

import std.algorithm;
import std.random;
import std.range;
import std.stdio;
import std.typecons;

class SpeciesClassificator {

    /// Constructor: Classify individuals into species.
    this( Phenotype[] individuals, float thresh) {
        this.thresh = thresh;
        // randomly pick first individual as prototype for species #0.
        prototypes ~= individuals[0];
        foreach(i; individuals) {
            auto best = bestMatch(i);
            if( best[0] > thresh ) {
                // individual i does not fit well into any species. Create new
                // species with i as prototype.
                prototypes ~= i.clone();
                i.species = cast(uint)prototypes.length - 1;
                writeln("Added species ", i.species);
            } else {
                // individual i fits into species
                i.species = best[1];
            }
        }
    }

    /// chose a random member of species as prototype for that species
    void updatePrototypes(Phenotype[] individuals) {
        foreach( idx, ref proto; prototypes ) {
            auto r = individuals.find!( i => i.species == idx ).array;
            proto = r[uniform(0,$)];
        }
    }

    /// shared fitness: Individuals' fitness is reduced by number of members in
    /// its species. This is done to prevent a single species becoming too
    /// dominant in numbers.
    float sharedFitness(Phenotype ind, Phenotype[] individuals, uint species) {
        // TODO: Save number of members in the species to prevent counting it again
        auto members = individuals.filter!(i=>i.species == species);
        return ind.fitness / members.walkLength();
    }

    /// get number of species
    uint numberOfSpecies() @property {
        return cast(uint)prototypes.length;
    }

private:

    /// Find prototype which has smallest distance to ind.
    /// Return distance and species index as tuple.
    auto bestMatch( Phenotype ind ) {
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

    Phenotype[] prototypes;
    float thresh;
}