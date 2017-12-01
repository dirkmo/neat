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
    Population pop = new Population(100, 2, 1, false);
    float[][] patterns = [
        [0, 0, 0],
        [0, 1, 1],
        [1, 0, 1],
        [1, 1, 0],
    ];

    foreach( gen; 0..100 ) {
        writeln("Generation ", gen);
        foreach( idx, ind; pop.individuals ) {
            write("\r  ind: ", idx);
            float totalError = 0;
            foreach( p; patterns ) {
                float output = ind.propagate( p[0..2] )[0];
                float error = p[2] - output;
                error = error * error;
                totalError += error;
            }
            ind.fitness = totalError;
        }
        writeln();
        pop.selection();
        pop.mutation();
        if( gen % 10 == 0 ) {
            pop.pool.resetMutationList();
        }
    }
}
