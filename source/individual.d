module neat.individual;

import neat.connection;
import neat.genepool;
import neat.node;
import neat.phenotype;

import std.algorithm;
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
        float[] newValues;
        newValues.length = pool.getNodeCount();
        newValues[] = 0.0f;
        // propagate
        foreach( c; cons ) {
            if( c.enabled ) {
                auto n1 = c.start;
                auto nid2 = c.end.id;
                newValues[nid2] += c.weight * sigmoid(n1.value);
            }
        }
        // copy new values into node.values
        foreach( n; nodes ) {
            n.value = newValues[n.id];
        }
        // return values of output nodes
        auto outputNodes = nodes[pool.inputs..(pool.inputs + pool.outputs)];
        float[] output;
        output.length = pool.outputs;
        foreach( i, on; outputNodes ) {
            output[i] = on.value;
        }
        return output;
    }

    const(float)[] propagate( float[] inputValues ) {                             
        assert( !pool.isRecurrent() );
        // copy values into input nodes
        foreach( i, inp; nodes[0..pool.inputs] ) {
            inp.value = inputValues[i];
        }
        foreach( n; nodes[pool.inputs..$] ) {
            n.value = 0.0f;
        }
        // propagate layer by layer
        foreach( lidx; 0..pool.getLayerCount()-1 ) {
            // get nodes in layer lidx
            auto layerNodes = nodes.filter!(n=>n.layerIndex==lidx)();
            foreach( n1; layerNodes ) {
                foreach( c; n1.getOutputConnections() ) {
                    if( c.enabled ) {
                        c.end.value += c.weight * sigmoid(n1.value);
                    }
                }
            }
        }
        // return values of output nodes
        auto outputNodes = nodes[pool.inputs..(pool.inputs + pool.outputs)];
        float[] output;
        output.length = pool.outputs;
        foreach( i, on; outputNodes ) {
            output[i] = on.value;
        }
        return output;
    }
    
private:
    float sigmoid( float x ) {
        return 1.0f / ( 1.0f + exp(-x) );
    }

}