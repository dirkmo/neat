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
    uint memberCount; // current member count
    uint nextGenMemberCount; // member count for next generation
    float fitness; // fitness sum over all members (no average)
    Individual prototype; // the species representative
}

class SpeciesClassificator {

    this(Individual[] individuals, float thresh) {
        this.thresh = thresh;
        update(individuals);
    }

    /// look for individuals without an assigned species
    /// and assign them to a species. If nearest species
    /// exceeds the threshold, a new species is created.
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

    /// pick new species prototypes and reassign individuals to
    /// nearest species.
    void reassign() {
        // choose new prototypes
        foreach( sp; species ) {
            sp.memberCount = 0;
            auto members = individuals.filter!(a => a.species == sp.index).array;
            sp.prototype = new Individual(members[uniform(0,$)]);
        }
        // assign members to species
        foreach(ind; individuals) {
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

    /// calculate fitness of species and total sum
    void calculateFitness() {
        uint count;
        totalFitness = 0;
        foreach(ref sp; species) {
            sp.fitness = 0;
            individuals.filter!(a => a.species == sp.index)
                       .each!(b => sp.fitness += b.fitness);
            totalFitness += sp.fitness;
            count += sp.memberCount;
        }
        assert(count == individuals.length, "ERROR: Not all individuals assigned a species!");
    }

    /// Calculate size of species regarding shared fitness for next generation.
    /// totalFitness has to be determined in prior by calculateFitness()
    void calculateNextGenSpeciesSize(uint popsize) {
        uint nextGenPopSize;
        foreach(sp; species) {
            sp.nextGenMemberCount = cast(uint)(popsize * sp.fitness / totalFitness);
            nextGenPopSize += sp.nextGenMemberCount;
        }
        // distribute remaining "free slots" over species sorted by fitness
        long rest = popsize - nextGenPopSize;
        if( rest > 0 ) {
            species.sort!( (a,b) => a.fitness < b.fitness );
            do {
                foreach(sp; species) {
                    sp.nextGenMemberCount++;
                    if( --rest == 0 ) {
                        break;
                    }
                }
            } while(rest > 0);
        }
        assert(rest == 0);
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
        assert( !s.empty );
        s.front.memberCount++;
    }

    Individual[] individuals;

    SpeciesData[] species;
    float thresh;

    float totalFitness;

    enum cExcess = 1.0f;
    enum cDisjunct = 1.0f;
    enum cWeight = 0.4f;
}