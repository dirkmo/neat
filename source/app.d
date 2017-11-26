import neat.genepool;
import neat.connection;
import neat.individual;
import neat.node;
import neat.population;

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

void printIndividual(Individual ind) {
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
        writefln("Con %s: %s -> %s ", c.innovation, c.start.id, c.end.id);
    }
}

void main()
{
    Population pop = new Population(2, 2, 1, false);
    printGenepool(pop.pool);
    writeln();
    Individual ind = pop.individuals[0];
    printIndividual(ind);
    writeln();

    
    Connection con = ind.cons[0];
    ind.mutateSplitUpConnection(con);
    printGenepool(pop.pool);
    writeln();
    printIndividual(pop.individuals[0]);
    writeln();

    Node n1, n2;
    n1 = ind.nodes[2];
    n2 = ind.nodes[3];
    ind.mutateAddConnection(n1, n2);
    printGenepool(pop.pool);
    writeln();
    printIndividual(pop.individuals[0]);
    writeln();

}
