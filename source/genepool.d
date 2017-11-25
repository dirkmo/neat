module neat.genepool;

import neat.connection;
import neat.node;

import std.algorithm;
import std.stdio;

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
        nodeAdded = false;
    }

    ///
    NodeGene createNewNodeGene() {
        auto newNode = new NodeGene(NodeGene.Type.hidden);
        nodeGenes ~= newNode;
        nodeAdded = true;
        return newNode;
    }

    ///
    bool mutateAddNewConGene( NodeGene n1, NodeGene n2, out ConGene newCon ) {
        writeln(__FUNCTION__);
        if( nodeAdded ) {
            updateNodesLayerIndex(null, 0);
        }
        
        if( !recurrent && n1.layerIndex >= n2.layerIndex ) {
            // connection would be recurrent
            writefln("Abort, connection would be recurrent." );
            return false;
        }

        if( n1.isInputToNode(n2, newCon) ) {
            // nodes are already connected.
            writefln("Connection gene %s already there", newCon);
            return true;
        }
        // create new connection gene
        newCon = new ConGene(n1, n2);
        writefln("Creating new congene %s", newCon);
        return true;
    }

    ///
    void mutateSplitUpConGene( ConGene oldCon, out ConGene con1, out ConGene con2 ) {
        writeln(__FUNCTION__);
        if( oldCon.innovation in splitUpConGenes ) {
            auto mutation = splitUpConGenes[oldCon.innovation];
            con1 = mutation.con1;
            con2 = mutation.con2;
            writefln("mutation already in genepool. Reusing con1: %s, con2: %s", con1, con2);
        } else {
            // split up connection gene
            NodeGene n1 = oldCon.start();
            NodeGene n2 = oldCon.end();
            // create new node
            NodeGene n3 = createNewNodeGene();
            // add connections for n1 --> n3 --> n2
            con1 = new ConGene(n1, n3);
            con2 = new ConGene(n3, n2);
            splitUpConGenes[oldCon.innovation] = SplitUpConGeneMutation(con1, con2);
        }
    }

    const(ConGene)[] getConGenes() const {
        return conGenes;
    }

private:
    uint inputs;
    uint outputs;
    bool recurrent;
    bool nodeAdded;

    struct SplitUpConGeneMutation {
        ConGene con1;
        ConGene con2;
    }

    SplitUpConGeneMutation[uint] splitUpConGenes; // con genes, that have been split up

    NodeGene[] nodeGenes;
    ConGene[] conGenes;
}