module neat;

import std.algorithm;
import std.container;
import std.random;
import std.range;
import std.stdio;

/+
Idee für Recurrent Networks:

Propagieren ist völlig Schichtunabhängig.
Es wird zwischen jeder vollständigen Propagierung der Zustand aller Nodes gespeichert.
ist an einem Node Ni ein rekurrentes Netz als Input angeschlossen, so wird der input des Netzes
aus der Aktivierung des Nodes Na der vorherigen Propagierung (Zeit t-1) verwendet.

Ni(t) += Na(t-1) * recurrent-net-weight
+/



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

    const(uint[]) getInputCons() { return inputCons; }
    const(uint[]) getOutputCons() { return outputCons; }

    void addOutputConnection( uint inno ) {
        if( !outputCons.canFind(inno) ) {
            outputCons ~= inno;
        }
    }

    void setLayer( int layer ) { layerIndex = layer; }
    int  getLayer() { return layerIndex; }

    Type getType() { return type; }

    uint getNodeId() { return nodeId; }

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

    /// perform split up mutation of a random connection gene
    void mutateSplitUpConGene(float probability) {
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
        if( uniform(0.0f, 1.0f) < probability ) {
            // choose two random nodes which are not connected
            // and connect them
            uint n1 = cast(uint)uniform(0, nodes.length);
            uint n2 = cast(uint)uniform(0, nodes.length);
            if( nodesAreConnected(n1, n2) ) {
                // already connected in phenotype
                return;
            }
            uint newCon;
            if( genepool.mutateAddNewConGene(n1, n2, newCon) ) {
                // get ConGene from gene pool. Will be added to gene pool
                // if it does not exist yet.
                const ConGene cg = genepool.getConGene(newCon);
                // create a phenotype fron ConGene
                ConPhenotype cp = new ConPhenotype(cg);
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
        // doesn't account for recurrent connections
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
    
    Genepool genepool;
    ConPhenotype[] cons;
    uint[] nodes;
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
        return newNode.nodeId;
    }

    bool mutateAddNewConGene(uint n1, uint n2, out uint newCon) {
        if( addedConGenes.canFind(newCon) ) {

            return true;
        }
        if( !recurrentAllowed && (nodeGenes[n1].getLayer() >= nodeGenes[n2].getLayer()) ) {
            // connection would be recurrent
            return false;
        }
        if( nodesAreConnected(n1, n2, newCon) ) {
            return false;
        }
        newCon = createNewConGene(n1, n2).getInnovation();
        addedConGenes ~= newCon; // save mutation
        return true;
    }

    void mutateSplitUpConGene( uint oldCon, ref uint con1, ref uint con2) {
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

    bool nodesAreConnected( uint n1, uint n2, ref uint con ) {
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
                // TODO: sanity check: es dürfen keine nodes mehr mit layer=-1 übrig sein.
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
        addedConGenes.length = 0;
    }

//  private:
    ConGene[] conGenes;
    NodeGene[] nodeGenes;
    uint innovationCount;
    bool recurrentAllowed;
    uint inputs;
    uint outputs;

    struct SplitUpConGeneMutation {
        uint con1;
        uint con2;
    }

    SplitUpConGeneMutation[uint] splitUpConGenes; // con genes, that have been split up
    uint[] addedConGenes; // con genes, that have been added
}
