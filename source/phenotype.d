module neat.phenotype;

import neat.connection;
import neat.genepool;
import neat.node;

import std.algorithm;

class Phenotype {
    this( Genepool pool, bool createConnections ) {
        this.pool = pool;
        if( createConnections ) {
            foreach( c; pool.getConGenes() ) {
                auto n1 = new Node(c.start());
                auto n2 = new Node(c.end());
                cons ~= new Connection(c, n1, n2);
                nodes ~= [n1, n2];
            }
        }
    }

    /// create an identical copy of this
    Phenotype clone() {
        Phenotype newp = new Phenotype(pool, false);
        foreach( c; cons ) {
            auto nid1 = c.input.id;
            auto nid2 = c.output.id;

            Node n1new, n2new;
            auto narr = nodes.find!(n=>n.id == nid1)();
            if ( narr.length == 0 ) {
                n1new = new Node(c.input.gene);
            } else {
                n1new = narr[0];
            }
            narr = nodes.find!(n=>n.id == nid2)();
            if( narr.length == 0 ) {
                n2new = new Node(c.output.gene);
            }
            auto cg = new Connection( c.gene, n1new, n2new );
            newp.nodes ~= [n1new, n2new];
            newp.cons ~= cg;
        }

        return newp;
    }

private:
    Genepool pool;
    Connection[] cons;
    Node[] nodes;
}