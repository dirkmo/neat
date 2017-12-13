module neat.species;

import neat.individual;
import neat.phenotype;

import std.algorithm;
import std.array;
import std.random;
import std.range;
import std.stdio;
import std.typecons;

struct SpeciesData {
    uint index;
    uint memberCount;
    uint nextMemberCount;
    float fitness;
    Individual prototype;
}

class SpeciesClassificator {

    this(Individual[] individuals, float thresh) {
        this.thresh = thresh;
        update(individuals);
    }

    void update(Individual[] individuals) {
        this.individuals = individuals;
        foreach(ind; individuals) {
            if( ind.species == uint.max ) {
                // no species assigned yet
                auto best = bestMatch(ind);
                if( best[0] < thresh ) {
                    ind.species = best[1];
                    addToSpecies(ind);
                } else {
                    // new species
                    species ~= SpeciesData(cast(uint)species.length, 1, 1, 0, new Individual(ind));
                }
            }
        }
    }

private:

    /// Find prototype which has smallest distance to ind.
    /// Return distance and species index as tuple.
    auto bestMatch( Individual ind ) {
        float bestDist = float.max;
        uint bestIdx = uint.max;
        foreach(idx, sp; species) {
            float dist = ind.distance(sp.prototype, cExcess, cDisjunct, cWeight );
            if( dist < bestDist ) {
                bestDist = dist;
                bestIdx = cast(uint)idx;
            }
        }
        writeln("Dist: ", bestDist);
        return tuple(bestDist, bestIdx);
    }

    void addToSpecies( Individual ind ) {
        auto s = species.find!( a => a.index == ind.species );
    }

    Individual[] individuals;

    SpeciesData[] species;
    float thresh;

    enum cExcess = 1.0f;
    enum cDisjunct = 1.0f;
    enum cWeight = 0.4f;
}