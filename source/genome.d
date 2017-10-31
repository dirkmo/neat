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
    @disable this();

    /// Create genome with inputs connected to all outputs,
    /// no hidden nodes
    this( uint inputs, uint outputs ) {
        this.inputs = inputs;
        this.outputs = outputs;
        foreach( i; 0..inputs ) {
            addNode( Node.Type.input );
        }
        foreach( o; 0..outputs ) {
            addNode( Node.Type.output );
        }
        foreach( i ; 0..inputs ) {
            foreach( o ; 0..outputs ) {
                addConnection( i, o );
            }
        }
    }

    // Clone genome g
    this( Genome g ) {
        inputs = g.inputs;
        outputs = g.outputs;

        // add nodes
        foreach( n ; g.nodes ) {
            addClonedNode(n);
        }
        foreach( c ; g.connections ) {
            addClonedConnection( c );
        }
    }

    Genome dup() {
        return new Genome(this);
    }

    /// add node with type, without connections
    Node addNode( Node.Type type ) {
        Node n = new Node(type);
        nodes[n.nodeId] = n;
        writeln("new node ", n.nodeId, " type: ", n.type);
        return n;
    }

    /// add new node with same type and id as n
    Node addClonedNode( Node n ) {
        Node newNode = new Node(n);
        nodes[n.nodeId] = newNode;
        return newNode;
    }

    /// add connection to map and update nodes connection list
    Connection addConnection( uint startNode, uint endNode ) {
        if( ! startNode in nodes ) {
            nodes[startNode] = new Node(Node.Type.hidden);
        }
        if( ! endNode in nodes ) {
            nodes[endNode] = new Node(Node.Type.hidden);
        }
        Connection c = new Connection( startNode, endNode );
        nodes[startNode].addOutput(c);
        nodes[endNode].addInput(c);
        connections[c.innovation] = c;
        return c;
    }

    Connection addClonedConnection( Connection c ) {
        Connection newCon = new Connection( c, c.inputNodeId, c.outputNodeId );
        nodes[c.inputNodeId].addOutput(newCon);
        nodes[c.outputNodeId].addInput(newCon);
        connections[c.innovation] = c;
        return c;
    }

    Genome crossOver( Genome g, bool gIsDominant ) {
        Genome offspring = new Genome( inputs, outputs );
        Connection[uint] p1 = connections.dup();
        Connection[uint] p2 = g.connections.dup();
        foreach( inno; p1.byKey ) {
            if( inno in p2 ) {
                // both parents have gene
                if( 0.5f < uniform(0.0f, 1.0f) ) {
                    Connection c = p1[inno];
                    offspring.addClonedNode( nodes[c.inputNodeId] );
                    offspring.addClonedNode( nodes[c.outputNodeId] );                    
                    offspring.addClonedConnection( c );
                } else {
                    Connection c = p2[inno];
                    offspring.addClonedNode( g.nodes[c.inputNodeId] );
                    offspring.addClonedNode( g.nodes[c.outputNodeId] );                    
                    offspring.addClonedConnection( c );              
                }
                p2.remove(inno);
            } else {
                // only p1 has gene. Take if p1 is dominant
                if( !gIsDominant ) {
                    Connection c = p1[inno];
                    offspring.addClonedNode( nodes[c.inputNodeId] );
                    offspring.addClonedNode( nodes[c.outputNodeId] );                    
                    offspring.addClonedConnection( c );
                }
            }
        }
        // now check remaining p2 genes if p2 is dominant
        if( gIsDominant ) {
            foreach( inno; p2.byKey ) {
                Connection c = p2[inno];
                offspring.addClonedNode( g.nodes[c.inputNodeId] );
                offspring.addClonedNode( g.nodes[c.outputNodeId] );                    
                offspring.addClonedConnection( c );
            }
        }
        return offspring;
    }

    Connection[uint] connections; // key is innovation number
    Node[uint] nodes; // key is node id
    uint inputs;
    uint outputs;
}