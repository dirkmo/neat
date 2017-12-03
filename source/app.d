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

    uint max = 1000;
    foreach( gen; 0..max ) {
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
                float error = abs(p[2] - output);
                totalError += error;
                //writefln("Output: %s", output);
            }
            ind.fitness = totalError + ind.cons.length / 20;
        }
        writeln();
        pop.selection();
        if(pop.first.fitness < 0.01 ) {
            break;
        }
        if( gen < max - 1 )
            pop.mutation();
    }
    float totalError = 0;
    foreach( pc, p; patterns ) {
        writeln("Pattern ", pc);
        float output = pop.individuals[0].propagate( p[0..2] )[0];
        writefln("Output: %s", output);
        float error = abs(p[2] - output);
        totalError += error;
        writefln("TotalError: %s", totalError);
    }
    
    printPhenotype(pop.first);

    for(int i=0; i<pop.individuals.length; i+=10) {
        auto ind = pop.individuals[i];
        writefln("%s: %s", i, ind.score());
    }
    
}
