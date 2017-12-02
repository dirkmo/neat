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
    Population pop = new Population(20, 2, 1, false);

    float[][] patterns = [
        [0, 0, 0],
        [0, 1, 1],
        [1, 0, 1],
        [1, 1, 0],
    ];

    foreach( gen; 0..10 ) {
        writeln("============================");
        writeln("Generation ", gen);
        writefln("Nodes: %s, Cons: %s, Layers: %s",
            pop.pool.nodeGenes.length, pop.pool.conGenes.length,
            pop.pool.getLayerCount());
        foreach( idx, ind; pop.individuals ) {
            //writefln("----------Ind: %s----------", idx);
            float totalError = 0;
            foreach( pc, p; patterns ) {
                //writeln("Pattern ", pc);
                float output = ind.propagate( p[0..2] )[0];
                float error = p[2] - output;                
                error = error * error;
                totalError += error;
                //writefln("Output: %s", output);
            }
            ind.fitness = totalError;
        }
        writeln();
        pop.selection();       
        if( gen < 48 ) 
            pop.mutation();
    }
    foreach( pc, p; patterns ) {
        writeln("Pattern ", pc);
        float output = pop.individuals[0].propagate( p[0..2] )[0];
        float error = p[2] - output;                
        error = error * error;
        writefln("Output: %s", output);
    }
    printPhenotype(pop.individuals[0]);

    ulong min, max;
    float mind = 1000, maxd = 0;
    foreach( idx, i; pop.individuals ) {
        if( idx == 0 ) continue;
        float d = pop.individuals[0].distance(i, 1, 1, 1);
        if( mind > d ) {
            min = idx; mind = d;
        }
        if( maxd < d ) {
            max = idx; maxd = d;
        }
        writefln("%s: %s", idx, d);
    }

}
