import neat.genepool;
import neat.connection;
import neat.individual;
import neat.node;
import neat.population;
import neat.phenotype;

import std.algorithm;
import std.container;
import std.math;
import std.random;
import std.range;
import std.stdio;


void printGenepool( Genepool pool ) {
    foreach(n; pool.nodeGenes ) {
        writefln("NodeGene %s (layer %s)", n.id, n.layerIndex);
        write("  Input ConGenes: ");
        foreach(c; n.getInputConGenes()) {
            writef("%s ", c.innovation);
        }
        writeln();
        write("  Output ConGenes: ");
        foreach(c; n.getOutputConGenes()) {
            writef("%s ", c.innovation);
        }
        writeln();
    }
    foreach(c; pool.getConGenes()) {
        writefln("ConGene %s: %s -> %s ", c.innovation, c.start.id, c.end.id);
    }
}

void printPhenotype(Phenotype ind) {
    foreach(n; ind.nodes ) {
        writefln("Node %s (layer %s)", n.id, n.layerIndex);
        write("  Input Cons: ");
        foreach(c; n.getInputConnections()) {
            writef("%s ", c.innovation);
        }
        writeln();
        write("  Output Cons: ");
        foreach(c; n.getOutputConnections()) {
            writef("%s ", c.innovation);
        }
        writeln();
    }
    foreach(c; ind.cons) {
        writefln("Con %s, %s -> %s, w: %s, %s",
                c.innovation, c.start.id, c.end.id, c.weight,
                c.enabled ? "enabled" : "disabled");
    }
}

void printState(Individual ind) {
    foreach(n; ind.nodes ) {
        writefln("Node %s: %s", n.id, n.value);
    }
    writeln();
}

void main()
{
    //string topology = "i0 i1 h0 b0 o0; i0-h0 i1-h0 i0-o0 i1-o0 b0-o0 b0-h0 h0-o0";
    string topology = "i0 i1 o0; i0-o0 i1-o0";

    Population pop = new Population(100, topology, false);

    float[][] patterns = [
        [0, 0, 0],
        [0, 1, 1],
        [1, 0, 1],
        [1, 1, 0],
    ];

    uint max = 600;
    enum FITNESS_TARGET = -0.05f;
    foreach( gen; 0..max ) {
        writeln("============================");
        writeln("Generation ", gen);
        const best = pop.best;
        writefln("Nodes: %s, Cons: %s, Layers: %s, best: %s (species %s)",
            pop.pool.nodeGenes.length, pop.pool.conGenes.length,
            pop.pool.getLayerCount(), best.fitness, best.species);
        foreach( idx, ind; pop.individuals ) {
            float totalError = 0;
            foreach( pc, p; patterns ) {
                float output = ind.propagate( p[0..2] )[0];
                const error = abs(output - p[2]);
                totalError += error;
            }
            ind.fitness = -totalError;
        }
        writeln();
        if(pop.best.fitness > FITNESS_TARGET ) {
            break;
        }
        pop.selection();
        pop.mutateWeights();
        if(gen % 5 == 0) {
            pop.mutateSplitUpConnections();
            pop.mutateAddConnections();
        }
    }

    float totalError = 0;

    auto best = pop.best;
    writefln("Best fitness: %s, species: %s", best.fitness, best.species);
    if( best.fitness > FITNESS_TARGET ) {
        printPhenotype(pop.best);
        foreach( pc, p; patterns ) {
            writef("Pattern %s: ", pc);
            float output = best.propagate( p[0..2] )[0];
            writefln("%s", output);
            const error = abs(output - p[2]);
            totalError += error;        
        }
        writefln("TotalError: %s", totalError);
    }
    
}
