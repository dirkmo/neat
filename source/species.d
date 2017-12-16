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
    float scale; // new species are scaled down at first
    Individual prototype; // the species representative

    float sharedFitness() {
        return fitness / memberCount;
    }
}

struct SpeciesRange {
    SpeciesData[] species;
    ulong idx;

    this(ref SpeciesData[] species) {
        this.species = species;
    }
    
    auto front() @property {
        return species[idx];
    }

    bool empty() @property {
        return idx >= species.length;
    }

    void popFront() {
        idx++;
    }

    SpeciesRange save() @property{
        return SpeciesRange(species);
    }
}

class SpeciesClassificator {

    this(Individual[] individuals, float thresh) {
        this.thresh = thresh;
        update(individuals);
        foreach(sp; species) {
            sp.scale = 1.0f;
        }
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
                    species ~= SpeciesData(
                        cast(uint)species.length, // index
                        1, // memberCount
                        1, // nextGenMemberCount
                        0.0, // fitness
                        NEW_SPECIES_SCALE, // scale
                        new Individual(ind) // prototype
                    );
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
                species ~= SpeciesData(cast(uint)species.length, 1, 1, 0.0, NEW_SPECIES_SCALE, new Individual(ind));
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
            sp.nextGenMemberCount = cast(uint)(popsize * sp.sharedFitness() * sp.scale / sharedFitness());
            if( sp.nextGenMemberCount < NEW_SPECIES_MEMBERS_MIN ) {
                sp.nextGenMemberCount = NEW_SPECIES_MEMBERS_MIN;
            } else if(sp.nextGenMemberCount > NEW_SPECIES_MEMBERS_MAX) {
                sp.nextGenMemberCount = NEW_SPECIES_MEMBERS_MAX;
            }
            nextGenPopSize += sp.nextGenMemberCount;
            sp.scale += sp.scale;
            if( sp.scale > 1.0f ) sp.scale = 1.0f;
        }
        // distribute remaining "free slots" over species sorted by fitness
        long rest = popsize - nextGenPopSize;
        if( rest > 0 ) {
            species.sort!( (a,b) => a.fitness > b.fitness );
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

    /// remove species, but not individuals
    void extinctSpecies(uint speciesIdx) {
        long idx = long.max;
        foreach(spidx, sp; species) {
            if(sp.index == speciesIdx) {
                idx = spidx;
                break;
            }
        }
        assert(idx != long.max);
    }

    SpeciesRange range() {
        return SpeciesRange(species);
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

    float sharedFitness() {
        return totalFitness / individuals.length;
    }

    Individual[] individuals;

    SpeciesData[] species;
    float thresh;

    float totalFitness;

    enum cExcess = 1.0f;
    enum cDisjunct = 1.0f;
    enum cWeight = 0.4f;
    enum NEW_SPECIES_MEMBERS_MIN = 1;
    enum NEW_SPECIES_MEMBERS_MAX = 10;
    enum NEW_SPECIES_SCALE = 0.1f;
}