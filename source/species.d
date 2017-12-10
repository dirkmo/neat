module neat.species;

import neat.phenotype;

import std.stdio;
import std.typecons;

class SpeciesClassificator {
    this() {
    }

    uint assignSpecies(Phenotype[] individuals, float thresh) {
        foreach(ref i; individuals) {
            auto best = bestMatch(i);
            if( best[0] > thresh ) {
                prototypes ~= i.clone();
                i.species = cast(uint)prototypes.length - 1;
                writeln("Added species ", i.species);
            } else {
                i.species = best[1];
            }
        }
        return cast(uint)prototypes.length;
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

    Phenotype[] prototypes;
}