module neat.genome;

import std.algorithm;
import std.random;
import std.stdio;

class Connection {
    this(uint inputNodeId, uint outputNodeId ) {
        this.inputNodeId = inputNodeId;
        this.outputNodeId = outputNodeId;
        weight = uniform(-1.0f, 1.0f);
        enabled = true;
        innovation = mutation.getNextInnovation();
    }

    this( Connection c, uint inputNodeId, uint outputNodeId ) {
        weight = c.weight;
        enabled = c.enabled;
        innovation = c.innovation;
        this.inputNodeId = inputNodeId;
        this.outputNodeId = outputNodeId;
    }
    
    void randomizeWeight() {
        weight = uniform(-1.0f, 1.0f);
    }
    
    void mutateWeights(float probability, float strength) {
        const auto r = uniform(0.0f, 1.0f);
        if( r < probability ) {
            if( r < probability/10 ) {
                weight = uniform(-1.0f, 1.0f);              
            } else {
                weight += weight * uniform(0.0f, strength);
            }
        }
    }

    uint inputNodeId; // connection gets its input from inputNodeId
    uint outputNodeId; // connection feeds data into outputNodeId
    float weight;
    bool enabled;
    uint innovation;
    
    static Mutation mutation;
}

class Node {
    enum Type {
        input,
        hidden,
        output
    }

    @disable this();

    this( Type type ) {
        this.type = type;
        nodeId = nodeCount++;
    }

    this( Node n ) {
        type = n.type;
        nodeId = n.nodeId;
    }

    this( Type type, uint nodeId ) {
        this.type = type;
        this.nodeId = nodeId;
    }

    bool isInput() { return type == Type.input; }
    bool isHidden() { return type == Type.hidden; }
    bool isOutput() { return type == Type.output; }

    void addInput( Connection c ) {
        input ~= c.innovation;
    }
    
    void addOutput( Connection c ) {
        output ~= c.innovation;
    }

    uint[] input; // innovation numbers of input nets
    uint[] output; // innovation numbers of output nets
    Type type;
    uint nodeId;
    static uint nodeCount;
}

class Mutation {
    static uint getNextInnovation() {
        return innovation++;
    }

    static void resetState() {

    }

    private static uint innovation;
}

class Genome {

    /// create empty genome
    this() {
    }

    /// Create genome with inputs connected to all outputs,
    /// no hidden nodes
    this( uint inputs, uint outputs ) {
        this.inputs = inputs;
        this.outputs = outputs;
        foreach( i; 0..inputs ) {
            addNodeWithId( Node.Type.input, i );
        }
        foreach( o; 0..outputs ) {
            addNodeWithId( Node.Type.output, o + inputs );
        }
        foreach( i ; 0..inputs ) {
            foreach( o ; 0..outputs ) {
                addConnection( i, inputs+o );
            }
        }
    }

    // Clone genome g
    this( Genome g ) {
        inputs = g.inputs;
        outputs = g.outputs;

        // add nodes
        foreach( n ; g.nodes ) {
            Node newNode = addClonedNode(n);
        }
        foreach( c ; g.connections ) {
            addClonedConnection( c );
        }
    }

    Genome dup() {
        return new Genome(this);
    }

    /// add node with type, new nodeId, without connections
    Node addNode( Node.Type type ) {
        Node n = new Node(type);
        nodes[n.nodeId] = n;
        return n;
    }

    Node addNodeWithId( Node.Type type, uint nodeId ) {
        Node n = new Node(type, nodeId);
        nodes[n.nodeId] = n;
        return n;
    }

    /// add new node with same type and id as n
    Node addClonedNode( Node n ) {
        if( ! (n.nodeId in nodes) ) {
            Node newNode = new Node(n);
            nodes[newNode.nodeId] = newNode;
            return newNode;
        }
        return nodes[n.nodeId];
    }

    /// add connection to map and update nodes connection list
    Connection addConnection( uint startNode, uint endNode ) {
        if( !(startNode in nodes) ) {
            nodes[startNode] = new Node(Node.Type.hidden, startNode);
        }
        if( !(endNode in nodes) ) {
            nodes[endNode] = new Node(Node.Type.hidden, endNode);
        }
        Connection c = new Connection( startNode, endNode );
        nodes[startNode].addOutput(c);
        nodes[endNode].addInput(c);
        connections[c.innovation] = c;
        return c;
    }

    Connection addClonedConnection( Connection c ) {
        Connection newCon = new Connection( c, c.inputNodeId, c.outputNodeId );
        nodes[newCon.inputNodeId].addOutput(newCon);
        nodes[newCon.outputNodeId].addInput(newCon);
        connections[newCon.innovation] = newCon;
        return newCon;
    }

    Genome crossOver( Genome g, bool gIsDominant ) {
        Genome offspring = new Genome();
        offspring.inputs = this.inputs;
        offspring.outputs = this.outputs;
        Connection[uint] p1 = connections.dup();
        Connection[uint] p2 = g.connections.dup();
        foreach( inno; p1.byKey ) {
            if( inno in p2 ) {
                // both parents have gene
                if( 0.5f < uniform(0.0f, 1.0f) ) {
                    Connection c = p1[inno];
                    Node n = offspring.addClonedNode( nodes[c.inputNodeId] );
                    n = offspring.addClonedNode( nodes[c.outputNodeId] );
                    Connection nc = offspring.addClonedConnection( c );
                } else {
                    Connection c = p2[inno];
                    Node n = offspring.addClonedNode( g.nodes[c.inputNodeId] );
                    n = offspring.addClonedNode( g.nodes[c.outputNodeId] );
                    Connection nc = offspring.addClonedConnection( c );              
                }
                p2.remove(inno);
            } else {
                // only p1 has gene. Take if p1 is dominant
                if( !gIsDominant ) {
                    Connection c = p1[inno];
                    offspring.addClonedNode( nodes[c.inputNodeId] );
                    offspring.addClonedNode( nodes[c.outputNodeId] );
                    Connection nc = offspring.addClonedConnection( c );
                }
            }

        }
        // now check remaining p2 genes if p2 is dominant
        if( gIsDominant ) {
            foreach( inno; p2.byKey ) {
                Connection c = p2[inno];
                offspring.addClonedNode( g.nodes[c.inputNodeId] );
                offspring.addClonedNode( g.nodes[c.outputNodeId] );                    
                Connection nc = offspring.addClonedConnection( c );
            }
        }
        return offspring;
    }

    Connection[uint] connections; // key is innovation number
    Node[uint] nodes; // key is node id
    uint inputs;
    uint outputs;
}