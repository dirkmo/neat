module neat.species;

import neat.individual;
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
        // pick first individual as prototype for species #0.
        prototypes ~= individuals[0].clone();
        foreach(i; individuals) {
            determineSpecies(i, individuals);
        }
    }

    ///
    this( Individual[] individuals, float thresh) {
        this(cast(Phenotype[])individuals, thresh);
    }

    /// chose a random member of species as prototype for that species
    void updatePrototypes(Phenotype[] individuals) {
        foreach( idx, ref proto; prototypes ) {
            auto r = individuals.find!( i => i.species == idx );
            proto = r[uniform(0,$)];
        }
    }

    /// pick new prototypes
    void updatePrototypes(Individual[] individuals) {
        updatePrototypes( cast(Phenotype[])individuals );
    }

    /// determine species for ind
    uint determineSpecies(Phenotype ind, Phenotype[] individuals) {
        auto best = bestMatch(ind);
        if( best[0] > thresh ) {
            // individual i does not fit well into any species. Create new
            // species with i as prototype.
            prototypes ~= ind.clone();
            ind.species = cast(uint)prototypes.length - 1;
            writeln("Added species ", ind.species);
        } else {
            // individual i fits into species
            ind.species = best[1];
        }
        return ind.species;
    }

    /// determine species for ind
    uint determineSpecies(Individual ind, Individual[] individuals) {
        return determineSpecies(cast(Phenotype)ind, cast(Phenotype[])individuals);
    }    

    /// shared fitness: Individuals' fitness is reduced by number of members in
    /// its species. This is done to prevent a single species becoming too
    /// dominant in numbers.
    float sharedFitness(Phenotype ind, Phenotype[] individuals, uint species) {
        // TODO: Save number of members in the species to prevent counting it again
        auto members = individuals.filter!(i=>i.species == species);
        return ind.fitness / members.walkLength();
    }

    /// get shared fitness for ind
    float sharedFitness(Individual ind, Individual[] individuals, uint species) {
        return sharedFitness(cast(Phenotype)ind, cast(Phenotype[])individuals, species);
    }

    /// calculate total fitness sum over all individuals
    float totalFitness(Individual[] individuals) {
        //writeln(__FUNCTION__);
        float totalFitness = 0.0f;
        foreach(ind; individuals) {
            totalFitness += sharedFitness(ind, individuals, ind.species);
        }
        //writeln("totalFitness: ", totalFitness);
        return totalFitness;
    }

    /// calculate fitness of species as a whole
    float speciesFitness(Individual[] individuals, uint species) {
        //writeln(__FUNCTION__);
        assert( species < prototypes.length );
        // calculate total fitness sum over all individuals
        float speciesFitness = 0.0f;
        individuals.filter!(i=>i.species == species)
                   .each!(m => speciesFitness += sharedFitness(m, individuals, species));
        //writeln("speciesFitness: ", speciesFitness);
        return speciesFitness;
    }

    /// delete species and reassign species ids
    /// this will delete only the prototype, but not the individuals belonging
    /// to that species.
    void extinction( uint[] species ) {
        foreach( sp; species ) {
            if( sp == 0 ) {
                prototypes = prototypes[1..$];
            } else if( sp == cast(uint)prototypes.length-1 ) {
                prototypes = prototypes[0..sp];
            } else {
                prototypes = prototypes[0..sp] ~ prototypes[sp+1..$];
            }
        }
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
            float dist = ind.distance(p, cExcess, cDisjunct, cWeight );
            if( dist < bestDist ) {
                bestDist = dist;
                bestIdx = cast(uint)idx;
            }
        }
        writeln("Dist: ", bestDist);
        return tuple(bestDist, bestIdx);
    }

    Phenotype[] prototypes;
    float thresh;

    enum cExcess = 1.0f;
    enum cDisjunct = 1.0f;
    enum cWeight = 0.4f;
}