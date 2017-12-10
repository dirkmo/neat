module neat.individual;

import neat.connection;
import neat.genepool;
import neat.node;
import neat.phenotype;

import std.algorithm;
import std.math;
import std.stdio;

///
class Individual : Phenotype {
    ///
    this( Genepool pool, bool createConnections ) {
        super( pool, createConnections );
    }

    ///
    const(float)[] propagateStep( in float[] inputValues ) {
        assert( inputValues.length == pool.inputs );
        // copy values into input nodes
        foreach( i, inp; nodes[0..pool.inputs] ) {
            inp.value = inputValues[i];
            //writefln("set node %s to %s", i, inp.value);
        }
        setBiasNodeValue(1.0f);        
        float[] newValues;
        newValues.length = pool.getNodeCount();
        newValues[pool.inputs..$] = 0.0f;
        newValues[0..pool.inputs] = inputValues;
        // propagate
        foreach( c; cons ) {
            if( c.enabled | true) {
                auto n1 = c.start;
                auto nid2 = c.end.id;
                newValues[nid2] += c.weight * activation(n1.value);
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
        foreach( i, n; nodes[0..pool.inputs] ) {
            n.value = inputValues[i];
        }
        foreach( n; nodes[pool.inputs..$] ) {
            n.value = 0.0f;
        }
        setBiasNodeValue(1.0f);
        //writeln("Layer count: ", pool.getLayerCount());
        // propagate layer by layer
        foreach( lidx; 0..pool.getLayerCount()-1 ) {
            // get nodes in layer lidx
            auto layerNodes = nodes.filter!(n=>n.layerIndex==lidx)();
            //writefln("Layer %s, Nodes: %s", lidx, layerNodes);
            foreach( n1; layerNodes ) {
                //writeln("Cons: ", n1.getOutputConnections() );
                foreach( c; n1.getOutputConnections() ) {
                    if( c.enabled ) {
                        float sig = activation(n1.value);
                        //writef("Node %s += con%s * Node %s: ", c.end.id, c.innovation, n1.id);
                        //writef("Node %s = %s + %s * %s = ", c.end.id, c.end.value, c.weight, n1.value);
                        c.end.value += c.weight * sig;
                        //writefln("%s", c.end.value);
                    }
                }
            }
        }
        // return values of output nodes
        auto outputNodes = nodes[pool.inputs..(pool.inputs + pool.outputs)];
        float[] output;
        output.length = pool.outputs;
        foreach( i, on; outputNodes ) {
            assert( on.type == NodeGene.Type.output, "Outputneuron ist kein output mehr." );
            output[i] = on.value;
        }
        return output;
    }


private:
    float activation( float x ) {
        //return 1.0f / ( 1.0f + exp(-x) );
        //return tanh(x);
        return x < 0 ? 0 : x;
        //return x;
    }

    void setBiasNodeValue( float val ) {
        // set bias node to value if it exists
        if( nodes.length > pool.inputs + pool.outputs &&
            nodes[pool.inputs+pool.outputs].type == NodeGene.Type.bias )
        {
              nodes[pool.inputs+pool.outputs].value = val;
        }
    }
}