import neat;
import std.math;
import std.random;
import std.stdio;

float sigmoid( float x ) {
    return 1.0f / ( 1.0f + exp(-x) );
}

void printGenes(Genepool gp) {
    foreach( cg; gp.conGenes ) {
        writefln( "congene %s: start: %s, end: %s", cg.getInnovation(), cg.getStartNodeId(), cg.getEndNodeId() );
    }
    foreach( ng; gp.nodeGenes ) {
        writefln("nodegene %s: type: %s, layer: %s", ng.getNodeId, ng.getType(), ng.getLayer());
        writef("  inputs(%s): ", ng.getInputCons().length);
        foreach( ic; ng.getInputCons() ) {
            writef("%s ", ic);
        }
        writeln();
        writef("  outputs (%s): ", ng.getOutputCons().length);
        foreach( oc; ng.getOutputCons() ) {
            writef("%s ", oc);
        }
        writeln();
    }
}

void main()
{
    Genepool genepool = new Genepool( 2, 1, true, false );
    uint old, c1, c2;
    genepool.mutateSplitUpConGene( old, c1, c2 );
    writefln("old: %s, c1: %s, c2: %s", old, c1, c2);
    printGenes( genepool );
}
