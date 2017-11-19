module neat;

import std.algorithm.comparison;
import std.algorithm;
import std.container;
import std.math;
import std.random;
import std.range;
import std.stdio;

class ConGene {
    @disable this();

    this( uint startId, uint endId, uint inno ) {
        startNodeId = startId;
        endNodeId = endId;
        innovation = inno;
    }

    void setStartNodeId( uint startNodeId ) { this.startNodeId = startNodeId; }
    uint getStartNodeId() const { return startNodeId; }

    void setEndNodeId( uint endNodeId ) { this.endNodeId = endNodeId; }    
    uint getEndNodeId() const { return endNodeId; }

    uint getInnovation() const { return innovation; }

  protected:
    uint startNodeId;
    uint endNodeId;
    uint innovation;
}

class ConPhenotype {
    @disable this();

    this( const ConGene gene ) {
        weight = uniform( -1.0f, 1.0f );
        this.gene = gene;
        enabled = true;
    }

    this( ConPhenotype cpt ) {
        this.gene = cpt.gene;
        this.weight = cpt.weight;
    }

    void mutateWeight( float probability, float strength ) {
        if( uniform(0.0f, 1.0f ) < probability ) {
            weight += strength * uniform(0.0f, 1.0f );
        }
    }

    const(ConGene) getConGene() { return gene; }

    float getWeight() { return weight; }

    bool enabled;

  private:
    const ConGene gene;
    float weight;
}

class NodeGene {
    enum Type { input, hidden, output }

    this( uint nodeId, Type type ) {
        this.nodeId = nodeId;
        this.type = type;
        if( type == Type.input ) {
            layerIndex = 0;
        }
    }

    void addInputConnection( uint inno ) {
        if( !inputCons.canFind(inno) ) {
            inputCons ~= inno;
        }
    }

    const(uint[]) getInputCons() const { return inputCons; }
    const(uint[]) getOutputCons() const { return outputCons; }

    void addOutputConnection( uint inno ) {
        if( !outputCons.canFind(inno) ) {
            outputCons ~= inno;
        }
    }

    void setLayer( int layer ) { layerIndex = layer; }
    int  getLayer() const { return layerIndex; }

    Type getType() const { return type; }

    uint getNodeId() const { return nodeId; }

  private:
    uint nodeId;
    Type type;
    int layerIndex = -1;
    uint[] inputCons;
    uint[] outputCons;
}

class Phenotype {

    this( Genepool pool, bool createConPhenotypes ) {
        genepool = pool;
        if( createConPhenotypes ){
            foreach( c; pool.conGenes ) {
                cons ~= new ConPhenotype( c );
            }
        }
        refreshNodesList();
    }

    /// perform split up mutation of a random connection phenotype
    void mutateSplitUpConGene(float probability) {
        writeln(__FUNCTION__);
        if( uniform(0.0f, 1.0f) < probability ) {
            uint oldCon = uniform(0, cast(uint)cons.length);
            if( !cons[oldCon].enabled ) {
                // gene is disabled, is already split up
                return;
            }
            uint con1inno, con2inno;
            uint oldConInno = cons[oldCon].getConGene().getInnovation();
            genepool.mutateSplitUpConGene( oldConInno, con1inno, con2inno);
            ConPhenotype ocp;
            assert(findConPhenotype( oldConInno, ocp ));
            // disable connection phenotype of old connection
            ocp.enabled = false;
            // create and add new connection phenotypes from genes
            ConPhenotype cp1 = new ConPhenotype( genepool.getConGene(con1inno) );
            ConPhenotype cp2 = new ConPhenotype( genepool.getConGene(con2inno) );
            cp1.weight = 1;
            cp2.weight = ocp.weight;
            cons ~= [ cp1, cp2 ];
            // add new node to nodes list
            nodes ~= cp1.getConGene().getEndNodeId();
            writefln("oldcon: %s, cp1: %s, cp2: %s", oldConInno, con1inno, con2inno);
        }
    }

    void mutateAddConPhenotype(float probability) {
        writeln(__FUNCTION__);
        if( uniform(0.0f, 1.0f) < probability ) {
            // choose two random nodes which are not connected
            // and connect them
            immutable n1 = nodes[uniform(0, nodes.length)];
            immutable n2 = nodes[uniform(0, nodes.length)];
            if( nodesAreConnected(n1, n2) ) {
                // already connected in phenotype
                writefln("nodes %s and %s already connected.", n1, n2);
                return;
            }
            uint newCon;
            if( genepool.mutateAddNewConGene(n1, n2, newCon) ) {
                // get ConGene from gene pool. Will be added to gene pool
                // if it does not exist yet.
                const ConGene cg = genepool.getConGene(newCon);
                // create a phenotype fron ConGene
                ConPhenotype cp = new ConPhenotype(cg);
                cons ~= cp;
            }
        }
    }

    /// mutate weights by probability with strength
    void mutateWeights(float probability, float strength) {
        // go over all connection phenotypes
        foreach(cp; cons) {
            cp.mutateWeight(probability, strength);
        }
    }

    /// create an identical copy of this
    Phenotype clone() {
        Phenotype newp = new Phenotype(genepool, false);
        foreach( c; cons ) {
            newp.cons ~= new ConPhenotype( c );
        }
        newp.nodes = nodes.dup();
        return newp;
    }

    Phenotype crossOver( Phenotype pt ) {
        Phenotype offspring = new Phenotype(genepool, false);
        auto p1 = SList!uint();
        auto p2 = SList!uint();
        
        foreach( c1; cons ) {
            p1.insertFront(c1.getConGene().getInnovation());
        }
        foreach( c2; pt.cons ) {
            p2.insertFront(c2.getConGene().getInnovation());
        }
        ConPhenotype cpt;
        foreach( i; p1 ) {
            auto range = p2[].find( i );
            if( !range.empty ) {
                // both parents have gene
                writefln("both: %s", i);
                if( uniform(0.0f, 1.0f) < 0.5f ) {
                    writeln("  take p1");
                    if( this.findConPhenotype(i, cpt) ) {
                        offspring.addConPhenotype(cpt);
                    }
                } else {
                    writeln("  take p2");
                    if( pt.findConPhenotype(i, cpt) ) {
                        offspring.addConPhenotype(cpt);
                    }
                }
                p2.linearRemove(range.take(1));
            } else {
                // only p1 has the gene
                writeln("p1: %i");
                if( this.findConPhenotype(i, cpt) ) {
                    offspring.addConPhenotype(cpt);
                }
            }
        }
        foreach( i; p2 ) {
            // only p2 has the gene            
            if( pt.findConPhenotype(i, cpt) ) {
                offspring.addConPhenotype(cpt);
            }
        }
        return offspring;
    }

    /// find phenotype of connection gene
    bool findConPhenotype( uint inno, ref ConPhenotype cpt ) {
        foreach(p; cons) {
            if( p.getConGene().getInnovation() == inno ) {
                cpt = p;
                return true;
            }
        }
        return false;
    }

    ConPhenotype addConPhenotype( ConPhenotype cpt ) {
        ConPhenotype newcpt = new ConPhenotype(cpt);
        cons ~= newcpt;
        return newcpt;
    }

    bool nodesAreConnected( uint n1, uint n2) {
        // doesn't account for backward connection
        foreach(cpt; cons) {
            const ConGene cg = cpt.getConGene();
            if( cg.getStartNodeId() == n1 && cg.getEndNodeId() == n2 ) {
                return true;
            }
        }
        return false;
    }

    void refreshNodesList() {
        nodes.length = 0;
        foreach( c; cons ) {
            uint n1 = c.getConGene().getStartNodeId();
            uint n2 = c.getConGene().getEndNodeId();
            if( !nodes.canFind(n1) ) {
                nodes ~= n1;
            }
            if( !nodes.canFind(n2) ) {
                nodes ~= n2;
            }
        }
    }

    uint getNodeCount() { return cast(uint)nodes.length; }

    //protected:
    Genepool genepool;
    ConPhenotype[] cons;
    uint[] nodes; // holds node gene ids 
}

class Genepool {
    @disable this();

    this( uint inputs, uint outputs, bool fullyConnected, bool enableRecurrentNets ) {
        this.inputs = inputs;
        this.outputs = outputs;
        // Input nodes have id 0..inputs-1
        foreach(i; 0..inputs) {
            nodeGenes ~= new NodeGene(cast(uint)nodeGenes.length, NodeGene.Type.input );
        }
        // Output nodes have id inputs..inputs+outputs-1
        foreach(i; 0..outputs) {
            nodeGenes ~= new NodeGene(cast(uint)nodeGenes.length, NodeGene.Type.output );
        }
        if( fullyConnected ) {
            foreach( i; 0..inputs ) {
                foreach( o; 0..outputs ) {
                    createNewConGene( getInputNodeId(i), getOutputNodeId(o) );
                }
            }
        }
        recurrentAllowed = enableRecurrentNets;
        // assign every node a layer index
        updateNodesLayerIndex(null, 0);
    }

    ConGene createNewConGene( uint inputNodeId, uint outputNodeId ) {
        ConGene cg = new ConGene( inputNodeId, outputNodeId, innovationCount++ );
        conGenes ~= cg;
        nodeGenes[inputNodeId].addOutputConnection(cg.getInnovation());
        nodeGenes[outputNodeId].addInputConnection(cg.getInnovation());
        return cg;
    }

    uint createNewHiddenNode() {
        auto newNode = new NodeGene(cast(uint)nodeGenes.length, NodeGene.Type.hidden );
        nodeGenes ~= newNode;
        nodeAdded = true;
        return newNode.nodeId;
    }

    /// add ConGene connection from n1 to n2 to genepool. newCon will hold the genes
    /// innovation number.
    /// if n1 and n2 are connected, the inno number of that connection will be returned.
    /// returns false, if connection would be recurrent and recurrent nets are disabled.
    bool mutateAddNewConGene(in uint n1, in uint n2, out uint newCon) {
        writeln(__FUNCTION__);
        if( nodeAdded ) {
            updateNodesLayerIndex(null, 0);
        }
        writefln("node %s (layer %s) to %s (layer %s):", n1, nodeGenes[n1].getLayer(),
                    n2, nodeGenes[n2].getLayer());

        if( !recurrentAllowed ) {
            if( nodeGenes[n1].getLayer() >= nodeGenes[n2].getLayer() ) {
                // connection would be recurrent
                writefln("Abort, connection would be recurrent." );
                return false;
            }
        }
        if( nodesAreConnected(n1, n2, newCon) ) {
            // nodes are already connected.
            writefln("Connection gene %s already there", newCon);
            return true;
        }
        // create new connection gene
        newCon = createNewConGene(n1, n2).getInnovation();
        writefln("Creating new congene %s", newCon);
        return true;
    }

    void mutateSplitUpConGene( in uint oldCon, out uint con1, out uint con2) {
        writeln(__FUNCTION__);
        if( oldCon in splitUpConGenes ) {
            // connection gene has already been split up
            auto mutation = splitUpConGenes[oldCon];
            con1 = mutation.con1;
            con2 = mutation.con2;
            writefln("mutation already in genepool. Reusing con1: %s, con2: %s", con1, con2);
        } else {
            // split up connection gene
            const uint n1 = conGenes[oldCon].startNodeId;
            const uint n2 = conGenes[oldCon].endNodeId;
            // create new node n3
            const uint n3 = createNewHiddenNode();
            // add connections for n1 --> n3 --> n2
            con1 = createNewConGene(n1, n3).getInnovation();
            con2 = createNewConGene(n3, n2).getInnovation();
            splitUpConGenes[oldCon] = SplitUpConGeneMutation(con1, con2); // save mutation
        }
    }

    bool nodesAreConnected( in uint n1, in uint n2, out uint con ) {
        auto ocons = nodeGenes[n1].outputCons;
        foreach( oc; ocons ) {
            if( oc == n2 ) {
                con = oc;
                return true;
            }
        }
        return false;
    }

    void updateNodesLayerIndex( NodeGene[] layerNodes, const int layerNumber ) {
        if( layerNumber == 0 ) {
            NodeGene[] list;
            // loop over all nodes
            foreach(n; nodeGenes) {
                // if node is input, set layer and add connected nodes to list
                if( n.getType == NodeGene.type.input ) {
                    n.setLayer(0);
                    // get every out connection of node
                    foreach(cp; n.getOutputCons()) {
                        // add nodes which are connected to this net
                        const uint nodeId = conGenes[cp].endNodeId;
                        NodeGene node = nodeGenes[nodeId];
                        // only add hidden nodes, skip output nodes, because they get the highest layer number
                        if( node.getType() == NodeGene.Type.hidden ) {
                            list ~= node;
                        }
                    }
                } else {
                    // mark layer as unknown for all nodes
                    n.setLayer(-1);
                }
            }
            // recursively call for next layer
            updateNodesLayerIndex( list, 1 );
        } else {
            if( layerNodes.length == 0 ) {
                // all hidden nodes are done.
                // now set layer of output nodes
                foreach(o; 0..outputs) {
                    nodeGenes[getOutputNodeId(o)].setLayer(layerNumber);
                }
                // recursion is done here
                // sanity check: there shall be no nodes left with layer=-1
                foreach(n; nodeGenes) {
                    assert( n.getLayer() > -1 );
                }
            } else {
                NodeGene[] list;
                foreach(n; layerNodes) {
                    n.setLayer(layerNumber);
                    // get every out connection of node
                    foreach(cp; n.getOutputCons()) {
                        // add nodes which are connected to this net, if they have no layer assigned yet and aren't output neurons
                        const uint nodeId = conGenes[cp].endNodeId;
                        NodeGene node = nodeGenes[nodeId];
                        if( (node.getType() == NodeGene.Type.hidden) && (node.getLayer() == -1) ) {
                            list ~= node;
                        }
                    }
                }
                // recursively call for next layer
                updateNodesLayerIndex( list, layerNumber + 1 );
            }
        }
        nodeAdded = false;
    }

    uint getInputNodeId( uint i ) {
        assert( i < inputs );
        return i;
    }

    uint getOutputNodeId( uint o ) {
        assert( o < outputs );
        return inputs + o;
    }

    const(ConGene) getConGene( uint innovation ) {
        assert( innovation < conGenes.length );
        return conGenes[innovation];
    }

    const(NodeGene) getNodeGene( uint nodeId ) {
        assert( nodeId < nodeGenes.length );
        return nodeGenes[nodeId];
    }

    void resetMutationMemory() {
        splitUpConGenes.clear();
    }

    bool isRecurrent() { return recurrentAllowed; }

    uint getNodeCount() { return cast(uint)nodeGenes.length; }
    uint getConGeneCount() { return cast(uint)conGenes.length; }

    uint getInputNodeCount() { return inputs; }
    uint getOutputNodeCount() { return outputs; }

    uint getLayerCount() {
        if( nodeAdded ) {
            updateNodesLayerIndex( null, 0 );
        }
        return getNodeGene(getOutputNodeId(0)).getLayer() + 1;
    }

//  private:
    ConGene[] conGenes;
    NodeGene[] nodeGenes;

    bool nodeAdded;

    uint innovationCount;
    bool recurrentAllowed;
    uint inputs;
    uint outputs;

    struct SplitUpConGeneMutation {
        uint con1;
        uint con2;
    }

    SplitUpConGeneMutation[uint] splitUpConGenes; // con genes, that have been split up
}

class Individual : Phenotype {
    this( Genepool pool, bool createConPhenotypes ) {
        super( pool, createConPhenotypes );
        nodeValues.length = pool.getNodeCount();
        nodeValues[] = 0.0f;
        updateConValues();
    }

    void updateConValues() {
        conValues.length = genepool.getConGeneCount();
        foreach( cpt; cons ) {
            conValues[cpt.getConGene().getInnovation()] = cpt.getWeight();
        }
    }

    float sigmoid( float x ) {
        return 1.0f / ( 1.0f + exp(-x) );
    }

    const(float)[] propagateStep( const float[] inputs ) {
        if( nodeValues.length !=  genepool.getNodeCount() ) {
            // genepool mutated. adding new nodes with value 0
            assert( nodeValues.length < genepool.getNodeCount() );
            float[] newNodes;
            newNodes.length = genepool.getNodeCount() - nodeValues.length;
            newNodes[] = 0.0f;
            nodeValues ~= newNodes;
        }
        float[] newValues;
        newValues.length = nodeValues.length;
        // copy input values
        const uint i1 = genepool.getInputNodeId( 0 );
        const uint i2 = genepool.getInputNodeId( genepool.getInputNodeCount() - 1 ) + 1;
        assert( i2-i1 == inputs.length );
        newValues[i1..i2] = inputs[];
        // set hidden and output neurons to 0.0f
        newValues[i2..$] = 0.0f;
        foreach(c; cons) {
            if( c.enabled ) {
                const(ConGene) cg = c.getConGene();
                const uint n1 = cg.getStartNodeId();
                const uint n2 = cg.getEndNodeId();
                const float w = c.getWeight();
                writefln("n1=%s, n2=%s", n1, n2 );
                writefln("n1: %s, n2: %s, w: %s, sig: %s", nodeValues[n1], nodeValues[n2], w, sigmoid(nodeValues[n1]) );
                newValues[n2] += sigmoid( nodeValues[n1] ) * w;
            }
        }
        nodeValues = newValues;
        // return output values
        uint o1 = genepool.getOutputNodeId( 0 );
        uint o2 = genepool.getOutputNodeId( genepool.getOutputNodeCount() - 1 ) + 1;
        return nodeValues[o1..o2];
    }

    const(float)[] propagate( float[] inputs ) {
        assert( !genepool.isRecurrent() );
        // TODO: dies in eigene Funktion auslagern. Muss auch nur gemacht
        // werden, wenn sich nodes[] geändert hat
        // setting up layer structure
        DList!(uint)[] layers;
        layers.length = genepool.getLayerCount();
        foreach(n; 0..genepool.getInputNodeCount()) {
            layers[0].insert(genepool.getInputNodeId(n));
        }
        foreach(c; cons) {
            auto cg = c.getConGene();
            uint nodeId = cg.getEndNodeId();
            uint l = genepool.getNodeGene(nodeId).getLayer();
            sortedInsert(layers[l], nodeId);
        }
        // now start propagating
        // copy input values
        const uint i1 = genepool.getInputNodeId( 0 );
        const uint i2 = genepool.getInputNodeId( genepool.getInputNodeCount() - 1 ) + 1;
        assert( i2-i1 == inputs.length );
        nodeValues[i1..i2] = inputs[];
        nodeValues[i2..$] = 0.0f;
        foreach( l; layers ) {
            // loop over all nodes in layer
            foreach( nid; l ) {
                // loop over all outgoing connections of node

            }

        }

        // return output values
        uint o1 = genepool.getOutputNodeId( 0 );
        uint o2 = genepool.getOutputNodeId( genepool.getOutputNodeCount() - 1 ) + 1;
        return nodeValues[o1..o2];
    }

    void sortedInsert( ref DList!(uint) list, uint nodeId ) {
        auto range = list[];
        while( !range.empty() ) {
            if( range.front == nodeId )
                return;
            if( range.front > nodeId ) {
                list.insertBefore( range, [nodeId] );
                return;
            }
            range.popFront();
        }
        list.insertBack( nodeId );
    }
    
    float[] nodeValues; // index is node id
    float[] conValues; // index is innovation number
}

class Population {
    this( uint popsize, uint inputs, uint outputs, bool recurrent ) {
        genepool = new Genepool(
            inputs, outputs,
            true /*fullyConnected*/,
            false /*enableRecurrentNets*/
        );
        pop.length = popsize;
        foreach( ref p; pop ) {
            p = new Individual( genepool, true /*createConPhenotype*/ );
        }
    }
    
    private:

    Genepool genepool;
    Individual[] pop;
}