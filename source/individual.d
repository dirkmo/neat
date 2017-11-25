module neat.individual;

import neat.connection;
import neat.genepool;
import neat.node;
import neat.phenotype;

class Individual : Phenotype {
    this( Genepool pool, bool createConPhenotypes ) {
        super( pool, createConPhenotypes );
    }
}