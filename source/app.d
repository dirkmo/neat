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
    Phenotype[] pop;
    foreach( i; 0..2 ) {
        pop ~= new Phenotype( genepool, true );
        pop[i].mutateWeights( 1.0, 1.0 );
        printPhenotype(pop[i]);
    }

    Phenotype child = pop[0].crossOver( pop[1] );
    printPhenotype(child);    
}
