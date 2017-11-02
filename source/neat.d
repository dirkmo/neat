module neat;

import std.algorithm;
import std.random;
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

    this( ConGene gene ) {
        weight = uniform( -1.0f, 1.0f );
        this.gene = gene;
        enabled = true;
    }

    void mutateWeight( float probability, float strength ) {
        if( uniform(0.0f, 1.0f ) < probability ) {
            weight += strength * uniform(0.0f, 1.0f );
        }
    }

    const(ConGene) getConGene() { return gene; }

    bool enabled;

  private:
    ConGene gene;
    float weight;
}

class NodeGene {
    enum Type { input, hidden, output }

    this( uint nodeId, Type type ) {
        this.nodeId = nodeId;
        this.type = type;
    }

    void addInputConnection( uint inno ) {
        inputCons ~= inno;
    }

    void addOutputConnection( uint inno ) {
        outputCons ~= inno;
    }

    Type getType() { return type; }

  private:
    uint nodeId;
    Type type;
    uint[] inputCons;
    uint[] outputCons;
}

class Phenotype {

    this( Genepool pool ) {
        genepool = pool;
    }

    void mutateSplitUpConGene(float probability) {
        if( uniform(0.0f, 1.0f) < probability ) {
            uint oldCon, con1, con2;
            genepool.mutateSplitUpConGene( oldCon, con1, con2);
            ConPhenotype ocp = findConPhenotype( oldCon );
            ConPhenotype cp1 = findConPhenotype( con1 );
            ConPhenotype cp2 = findConPhenotype( con2 );
            cp1.weight = 1, cp2.weight = ocp.weight;
            ocp.enabled = false;
        }
    }

    /// find phenotype of connection gene
    ConPhenotype findConPhenotype( uint inno ) {
        foreach(p; con) {
            if( p.getConGene().getInnovation() == inno ) {
                return p;
            }
        }
        return null;
    }

    Genepool genepool;
    ConPhenotype[] con;
}

class Genepool {
    @disable this();

    this( uint inputs, uint outputs, bool fullyConnected, bool enableConcurrentNets ) {
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
        concurrent = enableConcurrentNets;
    }

    ConGene createNewConGene( uint inputNodeId, uint outputNodeId ) {
        ConGene cg = new ConGene( inputNodeId, outputNodeId, innovationCount++ );
        conGenes ~= cg;
        return cg;
    }

    uint createNewHiddenNode() {
        auto newNode = new NodeGene(cast(uint)nodeGenes.length, NodeGene.Type.hidden );
        nodeGenes ~= newNode;
        return newNode.nodeId;
    }

    bool mutateAddNewConGene(ref uint newCon) {
        // choose two nodes which are not connected
        // and connect them
        uint retry;
        uint n1, n2;
        do {
            n1 = uniform( 0, cast(uint)nodeGenes.length );
            n2 = uniform( 0, cast(uint)nodeGenes.length );
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

  private:
    ConGene[] conGenes;
    NodeGene[] nodeGenes;
    uint innovationCount;
    bool concurrent;
    uint inputs;
    uint outputs;
}
