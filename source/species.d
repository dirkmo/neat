module neat.species;

import neat.individual;
import neat.phenotype;

import std.algorithm;
import std.array;
import std.random;
import std.range;
import std.stdio;
import std.string;
import std.typecons;

struct SpeciesData {
    uint index;
    uint memberCount; // current member count
    uint nextGenMemberCount; // member count for next generation
    float fitness; // fitness sum over all members (no average)
    uint age; // counts number of generations this species exists
    Individual prototype; // the species representative

    static uint globalSpeciesCounter;

    enum AGE_MATURE = 10;


    this( Individual prototype ) {
        index = globalSpeciesCounter++;
        memberCount = 1;
        nextGenMemberCount = 1;
        fitness = 0;
        this.prototype = prototype;
    }
    
    /// fitness adjusted for species size and age
    @property float sharedFitness() {
        const scale = cast(float)age / AGE_MATURE;
        return fitness * scale / memberCount;
    }

    @property bool isNew() { return age < AGE_MATURE; }
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

    this(Individual[] individuals, float thresh, uint speciesCountMax) {
        this.thresh = thresh;
        this.speciesCountMax = speciesCountMax;
        update(individuals);
        foreach( ref sp; species ) {
            sp.age = SpeciesData.AGE_MATURE;
        }
    }

    /// look for individuals without an assigned species
    /// and assign them to a species. If nearest species
    /// exceeds the threshold, a new species is created.
    void update(Individual[] individuals) {
        writeln(__FUNCTION__);
        this.individuals = individuals;
        foreach(ind; individuals) {
            if( ind.species == uint.max ) {
                assignIndividual( ind );
            }
        }
        countSpeciesMembers();        
    }

    /// pick new species prototypes
    void pickNewPrototypes(Individual[] individuals) {
        writeln(__FUNCTION__);
        this.individuals = individuals;
        // choose new prototypes
        foreach( sp; species ) {
            sp.memberCount = 0;
            auto members = individuals.filter!(a => a.species == sp.index).array;
            if( members.length ) {
                sp.prototype = new Individual(members[uniform(0,$)]);
            } else {
                extinctSpecies(sp.index);
            }
        }
        // assign members to a species
//      foreach(ind; individuals) {
//          assignIndividual(ind);
//      }
        countSpeciesMembers();
    }

    /// calculate fitness of species and total sum
    void calculateFitness() {
        writeln(__FUNCTION__);
        uint count;
        totalFitness = 0;
        foreach(ref sp; species) {
            sp.fitness = 0;
            individuals.filter!(a => a.species == sp.index)
                       .each!(b => sp.fitness += b.fitness);
            totalFitness += sp.sharedFitness;
            count += sp.memberCount;
        }
        writefln("Count: %s, length: %s", count, individuals.length);
        writefln("totalFitness: %s", sharedFitness());
        assert(count == individuals.length, "ERROR: Not all individuals assigned a species!");
    }

    /// Calculate size of species regarding shared fitness for next generation.
    /// totalFitness has to be determined in prior by calculateFitness()
    void calculateNextGenSpeciesSize(uint popsize) {
        writeln(__FUNCTION__);
        uint nextGenPopSize;
        foreach(ref sp; species) {
            writefln("Species %s: MemberCount: %s, age: %s", sp.index, sp.memberCount, sp.age);
            writefln("  sp.sharedFitness: %s, sharedFitness: %s", sp.sharedFitness(), sharedFitness());
            sp.nextGenMemberCount = cast(uint)(popsize * sp.sharedFitness / sharedFitness());
            if( sp.isNew ) {
                if( sp.nextGenMemberCount < 1 ) sp.nextGenMemberCount = 1;
            }
            writefln("  nextGenMemberCount: %s", sp.nextGenMemberCount);
            nextGenPopSize += sp.nextGenMemberCount;
        }
        long rest = cast(long)popsize - nextGenPopSize;
        if( rest > 0 ) {
            writefln("rest: %s, popsize: %s, nextGenPopSize: %s", rest, popsize, nextGenPopSize);
            // distribute remaining "free slots" over species sorted by fitness
            species.sort!( (a,b) => a.fitness > b.fitness );
            do {
                foreach(ref sp; species) {
                    sp.nextGenMemberCount++;
                    if( --rest <= 0 ) {
                        break;
                    }
                }
            } while(rest > 0);
        } else if( rest < 0 ) {
            // remove excess members
            species.sort!( (a,b) => a.nextGenMemberCount > b.nextGenMemberCount );
            while( rest < 0 ) {
                foreach(ref sp; species) {
                    if( rest < 0 && !sp.isNew ) {
                        rest++;
                        sp.nextGenMemberCount--;
                    }
                }
            }
        }
        foreach(sp; species) {
            writefln("  nextGenMemberCount: %s", sp.nextGenMemberCount);
        }
        writefln("rest: %s, popsize: %s, nextGenPopSize: %s", rest, popsize, nextGenPopSize);
        assert(rest == 0);
    }

    /// remove species, but not individuals
    void extinctSpecies(uint speciesIdx) {
        foreach(spidx, sp; species) {
            if(sp.index == speciesIdx) {
                writefln("Species %s extinct");
                species.remove(spidx);
                return;
            }
        }
        throw new Exception(format("Species %s not found", speciesIdx));
    }

    ///
    SpeciesRange range() {
        return SpeciesRange(species);
    }

    ///
    @property uint numberOfSpecies() const {
        return cast(uint)species.length;
    }

    /// get species by id (param species is not an array index)
    void getSpeciesData(uint species, out SpeciesData speciesData) {
        foreach( ref sp; this.species ) {
            if( sp.index == species ) {
                speciesData = sp;
            }
        }
        throw new Exception(format("Species %s does not exist", species));
    }

    void age() {
        foreach(ref sp; species)
            sp.age++;
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
        return tuple(bestDist, bestIdx);
    }

    void countSpeciesMembers() {
        writeln(__FUNCTION__);
        uint[uint] count;
        foreach( ind; individuals ) {
            if( ind.species == uint.max ) {
                assignIndividual(ind);
            }
            if( ind.species in count ) {
                count[ind.species]++;
            } else {
                count[ind.species] = 1;
            }
        }
        writeln("Vorher: ", count);
        foreach( ref sp; species ) {
            writef("Species %s ", sp.index);
            stdout.flush();
            sp.memberCount = count[sp.index];
            writefln("has %s members", sp.memberCount);
            count.remove( sp.index );
        }
        assert( count.length == 0 );
    }

    void assignIndividual( Individual ind ) {
        auto best = bestMatch(ind);
        if( (best[0] < thresh) || (species.length == speciesCountMax) ) {
            ind.species = best[1];
        } else {
            // new species
            species ~= SpeciesData(new Individual(ind));
            ind.species = species[$-1].index;
            writefln("New species %s created.", species[$-1].index);
        }
    }

    float sharedFitness() {
        return totalFitness/* / individuals.length*/;
    }

    Individual[] individuals;

    SpeciesData[] species;
    float thresh;

    float totalFitness;

    uint speciesCountMax;

    enum cExcess = 1.0f;
    enum cDisjunct = 1.0f;
    enum cWeight = 0.4f;
}