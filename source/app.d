import neat;
import std.math;
import std.random;
import std.stdio;

void printGenes(Genepool gp) {
    foreach( cg; gp.conGenes ) {
        writefln( "congene %s: start: %s, end: %s", cg.getInnovation(), cg.getStartNodeId(), cg.getEndNodeId() );
    }
    foreach( ng; gp.nodeGenes ) {
        writefln("nodegene %s: type: %s, layer: %s", ng.getNodeId, ng.getType(), ng.getLayer());
        writef("  input cons(%s): ", ng.getInputCons().length);
        foreach( ic; ng.getInputCons() ) {
            writef("%s ", ic);
        }
        writeln();
        writef("  output cons (%s): ", ng.getOutputCons().length);
        foreach( oc; ng.getOutputCons() ) {
            writef("%s ", oc);
        }
        writeln();
    }
}

void printPhenotype( Phenotype p ) {
    foreach( c; p.cons ) {
        writefln( "congene %s: start: %s, end: %s, w: %s, %s", c.getConGene().getInnovation(),
                c.getConGene().getStartNodeId(), c.getConGene().getEndNodeId(), c.getWeight,
                c.enabled ? "enabled" : "disabled" );
    }
}

void main()
{
    Genepool genepool = new Genepool( 2, 1, true, false );
    Individual[] ind;
    foreach( i; 0..1 ) {
        writeln("pt ", i);
        ind ~= new Individual( genepool, true );
        ind[i].mutateWeights( 1.0, 1.0 );
        foreach(abc ; 0..3) {
            ind[i].mutateSplitUpConGene( 1.0f );
        }
        printPhenotype(ind[i]);
    }

    float[] input = [0.0f, 0.0f];
    auto output = ind[0].propagateStep( input );
    foreach(o; output) {
        writefln("Output: %s", o);
    }
}
