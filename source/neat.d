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
    }

    /// perform split up mutation of a random connection gene
    void mutateSplitUpConGene(float probability) {
        if( uniform(0.0f, 1.0f) < probability ) {
            uint oldCon, con1, con2;
            genepool.mutateSplitUpConGene( oldCon, con1, con2);
            ConPhenotype ocp;
            assert(findConPhenotype( oldCon, ocp ));
            ConPhenotype cp1 = new ConPhenotype( genepool.getConGene(con1) );
            ConPhenotype cp2 = new ConPhenotype( genepool.getConGene(con2) );
            ocp.enabled = false;
            cp1.weight = 1;
            cp2.weight = ocp.weight;
            cons ~= [ cp1, cp2 ];
        }
    }

    void mutateAddConPhenotype(float probability) {
        if( uniform(0.0f, 1.0f) < probability ) {
            uint newCon;
            if( genepool.mutateAddNewConGene(newCon) ) {
                // new ConGene added to gene pool
                // create a phenotype
                const ConGene cg = genepool.getConGene(newCon);
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

    Genepool genepool;
    ConPhenotype[] cons;
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

    bool mutateAddNewConGene(ref uint newCon) {
        // choose two random nodes which are not connected
        // and connect them
        // TODO: evtl. besser nur ein Versuch
        uint retry;
        uint n1, n2;
        do {
            n1 = uniform( 0, cast(uint)nodeGenes.length );
            n2 = uniform( 0, cast(uint)nodeGenes.length );
            writefln("%s %s", n1, n2);
            if( !recurrentAllowed && (nodeGenes[n1].getLayer() >= nodeGenes[n2].getLayer()) ) {
                // connection would be recurrent
                continue;
            }
            retry++;
        } while(nodesAreConnected( n1, n2 ) && (retry < 100));
        if( retry == 100 ) {
            return false;
        }
        newCon = createNewConGene(n1, n2).getInnovation();
        return true;
    }

    void mutateSplitUpConGene( ref uint oldCon, ref uint con1, ref uint con2) {
        oldCon = uniform(0, cast(uint)conGenes.length);
        const uint n1 = conGenes[oldCon].startNodeId;
        const uint n2 = conGenes[oldCon].endNodeId;
        // create new node n3
        const uint n3 = createNewHiddenNode();
        // add connections for n1 --> n3 --> n2
        con1 = createNewConGene(n1, n3).getInnovation();
        con2 = createNewConGene(n3, n2).getInnovation();
    }

    bool nodesAreConnected( uint n1, uint n2 ) {
        auto ocons = nodeGenes[n1].outputCons;
        foreach( oc; ocons ) {
            if( oc == n2 ) {
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

//  private:
    ConGene[] conGenes;
    NodeGene[] nodeGenes;
    uint innovationCount;
    bool recurrentAllowed;
    uint inputs;
    uint outputs;
}
