module neat.individual;

import neat.connection;
import neat.genepool;
import neat.node;
import neat.phenotype;

import std.math;

///
class Individual : Phenotype {
    ///
    this( Genepool pool, bool createConnections ) {
        super( pool, createConnections );
    }

    ///
    const(float)[] propagateStep( in float[] inputValues ) {
        // copy values into input nodes
        foreach( i, inp; nodes[0..pool.inputs] ) {
            inp.value = inputValues[i];
        }


        float[] outputs;
        outputs.length = pool.outputs;
        return outputs;
    }
    
private:
    float sigmoid( float x ) {
        return 1.0f / ( 1.0f + exp(-x) );
    }

}