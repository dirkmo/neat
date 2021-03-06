module neat.genepool;

import neat.connection;
import neat.node;

import std.algorithm;
import std.conv;
import std.stdio;
import std.string;

///
class Genepool {
    ///
    this( uint _inputs, uint _outputs, bool addBias, bool fullyConnected, bool enableRecurrentNets ) {
        this._inputs = _inputs;
        this._outputs = _outputs;
        // Input nodes have id 0.._inputs-1
        foreach(i; 0.._inputs) {
            nodeGenes ~= new NodeGene( NodeGene.Type.input );
        }
        // Output nodes have id _inputs.._inputs+_outputs-1
        foreach(i; 0.._outputs) {
            nodeGenes ~= new NodeGene( NodeGene.Type.output );
        }
        if( addBias ) {
            createBiasNodeGene();
            // bias neuron has id inputs + outputs
        }
        if( fullyConnected ) {
            auto inputNodes = nodeGenes.filter!(n=>n.type == NodeGene.Type.input)();
            auto outputNodes = nodeGenes.filter!(n=>n.type == NodeGene.Type.output)();
            foreach( i; inputNodes ) {
                foreach( o; outputNodes ) {
                    auto con = new ConGene(i, o);
                    conGenes ~= con;
                }
            }
        }
        recurrent = enableRecurrentNets;
        // assign every node a layer index
        updateNodesLayerIndex(null, 0);
    }

    /// topology: string that defines start topology.
    /// example: i0 i1 h1 o0 b0 ; i0-h1 i1-h1 h1-o1 b0-h1
    /// i: input, h: hidden, b: bias, o: output
    this( string topology, bool enableRecurrentNets ) {
        // first: adding nodes in correct order
        assert(topology.canFind(';'));
        auto split = topology.until(';').to!string().splitter(' ');
        auto input = split.filter!(s=>s[0] == 'i')();
        auto hidden = split.filter!(s=>s[0] == 'h')();
        auto output = split.filter!(s=>s[0] == 'o')();
        auto bias = split.filter!(s=>s[0] == 'b')();
        NodeGene[string] nodeMap;
        foreach(n; input) {
            auto node = new NodeGene( NodeGene.Type.input );
            nodeGenes ~= node;
            nodeMap[n] = node;
            _inputs++;
        }
        foreach(n; output) {
            auto node = new NodeGene( NodeGene.Type.output );
            nodeGenes ~= node;
            nodeMap[n] = node;
            _outputs++;
        }
        foreach(n; bias) {
            auto node = new NodeGene( NodeGene.Type.bias );
            nodeGenes ~= node;
            nodeMap[n] = node;
        }
        foreach(n; hidden) {
            auto node = new NodeGene( NodeGene.Type.hidden );
            nodeGenes ~= node;
            nodeMap[n] = node;
        }

        // create connections
        split = topology.find(';').to!string().splitter(' ');
        split.popFront();
        foreach( c; split ) {
            assert(c[2]=='-');
            NodeGene n1 = nodeMap[c[0..2]];
            NodeGene n2 = nodeMap[c[3..$]];
            auto con = new ConGene(n1, n2);
            conGenes ~= con;
        }

        recurrent = enableRecurrentNets;
        updateNodesLayerIndex(null, 0);
    }

    ///
    void updateNodesLayerIndex( NodeGene[] layerNodes, const int layerIndex ) {
        //writeln("LayerIndex: ", layerIndex);
        if( layerIndex == 0 ) {
            NodeGene[] list;
            // loop over all nodes
            foreach( n; nodeGenes ) {
                // if node is input, set layer and add connected nodes to list
                if( n.type == NodeGene.Type.input || n.type == NodeGene.Type.bias ) {
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
                layerCount = layerIndex + 1;        
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

    NodeGene createBiasNodeGene() {
        auto newNode = new NodeGene(NodeGene.Type.bias);
        nodeGenes ~= newNode;
        nodeAdded = true;
        return newNode;
    }

    ///
    bool mutateAddNewConGene( NodeGene n1, NodeGene n2, out ConGene newCon ) {
        //writeln(__FUNCTION__);
        if( n2.type == NodeGene.Type.bias ) {
            // bias neurons cannot have inputs, only outputs.
            return false;
        }

        if( nodeAdded ) {
            updateNodesLayerIndex(null, 0);
        }
        if( !recurrent && n1.layerIndex >= n2.layerIndex ) {
            // neural net is not recurrent, but connection would be recurrent
            //writefln("Abort, connection would be recurrent." );
            return false;
        }

        if( n1.isInputToNode(n2, newCon) ) {
            // nodes are already connected.
            //writefln("Connection gene %s already there", newCon);
            return true;
        }
        // create new connection gene
        newCon = new ConGene(n1, n2);
        conGenes ~= newCon;
        //writefln("Creating new congene %s", newCon);
        return true;
    }

    ///
    void mutateSplitUpConGene( ConGene oldCon, out ConGene con1, out ConGene con2 ) {
        //writeln(__FUNCTION__);
        if( oldCon.innovation in splitUpConGenes ) {
            auto mutation = splitUpConGenes[oldCon.innovation];
            con1 = mutation.con1;
            con2 = mutation.con2;
            //writefln("mutation reuse for %s: cg1: %s, cg2: %s", oldCon, oldCon.innovation, con1, con2);
        } else {
            // split up connection gene
            NodeGene n1 = oldCon.start();
            NodeGene n2 = oldCon.end();
            // create new node
            NodeGene n3 = createNewNodeGene();
            // add connections for n1 --> n3 --> n2
            con1 = new ConGene(n1, n3);
            con2 = new ConGene(n3, n2);
            conGenes ~= [con1, con2];
            splitUpConGenes[oldCon.innovation] = SplitUpConGeneMutation(con1, con2);
            //writefln("new connection for %s: cg1: %s, cg2: %s", oldCon, con1, con2);
        }
    }

    void resetMutationList() {
        splitUpConGenes.clear();
    }

    ConGene[] getConGenes() { return conGenes; }
    uint getNodeCount() const { return cast(uint)nodeGenes.length; }

    uint inputs() const @property { return _inputs; }
    uint outputs() const @property { return _outputs; }
    bool isRecurrent() const { return recurrent; }
    uint getLayerCount() const { return layerCount; }

//private:
    uint _inputs;
    uint _outputs;
    bool recurrent;
    uint layerCount;
    bool nodeAdded;

    struct SplitUpConGeneMutation {
        ConGene con1;
        ConGene con2;
    }

    SplitUpConGeneMutation[uint] splitUpConGenes; // con genes, that have been split up

    NodeGene[] nodeGenes;
    ConGene[] conGenes;
}