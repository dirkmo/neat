module neat.genepool;

import neat.connection;
import neat.node;

import std.algorithm;

///
class Genepool {
    ///
    this( uint inputs, uint outputs, bool fullyConnected, bool enableRecurrentNets ) {
        this.inputs = inputs;
        this.outputs = outputs;
        // Input nodes have id 0..inputs-1
        foreach(i; 0..inputs) {
            nodeGenes ~= new NodeGene( NodeGene.Type.input );
        }
        // Output nodes have id inputs..inputs+outputs-1
        foreach(i; 0..outputs) {
            nodeGenes ~= new NodeGene( NodeGene.Type.output );
        }
        if( fullyConnected ) {
            auto inputNodes = nodeGenes.filter!(n=>n.type == NodeGene.Type.input)();
            auto outputNodes = nodeGenes.filter!(n=>n.type == NodeGene.Type.output)();
            foreach( i; inputNodes ) {
                foreach( o; outputNodes ) {
                    conGenes ~= new ConGene(i, o);
                }
            }
        }
        recurrent = enableRecurrentNets;
        // assign every node a layer index
        updateNodesLayerIndex(null, 0);
    }

    ///
    void updateNodesLayerIndex( NodeGene[] layerNodes, const int layerIndex ) {
        if( layerIndex == 0 ) {
            NodeGene[] list;
            // loop over all nodes
            foreach( n; nodeGenes ) {
                // if node is input, set layer and add connected nodes to list
                if( n.type == NodeGene.Type.input ) {
                    n.layerIndex = 0;
                    // get every out connection of node
                    foreach( c; n.getOutputConGenes() ) {
                        // add nodes which are connected to this net
                        if( c.end().type == NodeGene.Type.hidden ) {
                            list ~= c.end();
                        }
                    }
                } else {
                    // mark layer as unknown
                    n.layerIndex = -1;
                }
            }
            updateNodesLayerIndex( list, 1 );
        } else {
            if( layerNodes.length == 0 ) {
                // all hidden nodes are done.
                // now set layer of output nodes
                auto outputNodes = nodeGenes.filter!(n=>n.type==NodeGene.Type.output)();
                foreach( n; outputNodes ) {
                    n.layerIndex = layerIndex;
                }
                // recursion is done here
                // sanity check: there shall be no nodes left with layer=-1
                foreach(n; nodeGenes) {
                    assert( n.layerIndex > -1 );
                }                
            } else {
                NodeGene[] list;
                foreach( n; layerNodes ) {
                    n.layerIndex = layerIndex;
                    // get every out connection of node
                    foreach(c; n.getOutputConGenes()) {
                        // add nodes which are connected to this net, if they have no layer assigned yet and aren't output neurons
                        NodeGene ng = c.end();
                        if( (ng.type == NodeGene.Type.hidden) && (ng.layerIndex == -1) ) {
                            list ~= ng;
                        }
                    }
                }
                // recursively call for next layer
                updateNodesLayerIndex( list, layerIndex + 1 );
            }
        }
    }

private:
    uint inputs;
    uint outputs;
    bool recurrent;

    NodeGene[] nodeGenes;
    ConGene[] conGenes;
}