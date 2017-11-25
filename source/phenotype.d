module neat.phenotype;

import neat.connection;
import neat.genepool;
import neat.node;

class Phenotype {
    this( Genepool pool, bool createConnections ) {
        this.pool = pool;
        if( createConnections ) {
            foreach( c; pool.getConGenes() ) {
                auto n1 = new Node(c.start());
                auto n2 = new Node(c.end());
                cons ~= new Connection(c, n1, n2);
            }
        }
    }

private:
    Genepool pool;
}