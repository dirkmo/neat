module neat.genome;

import std.algorithm;
import std.random;
import std.stdio;

class Connection {
    this(uint inputNode, uint outputNode ) {
        this.inputNode = inputNode;
        this.outputNode = outputNode;
        weight = uniform(-1.0f, 1.0f);
        enabled = true;
        innovation = mutation.getNextInnovation();
    }

    this( Connection c, uint inputNode, uint outputNode ) {
        weight = c.weight;
        enabled = c.enabled;
        innovation = c.innovation;
        this.inputNode = inputNode;
        this.outputNode = outputNode;
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

    uint inputNode; // connection gets its input from inputNode
    uint outputNode; // connection feeds data into outputNode
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

    bool isInput() { return type == Type.input; }
    bool isHidden() { return type == Type.hidden; }
    bool isOutput() { return type == Type.output; }

    void addInput( Connection c ) {
        input ~= c.innovation;
    }
    
    void addOutput( Connection c ) {
        output ~= c.innovation;
    }

    uint[] input;
    uint[] output;
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
    @disable this();
    this( uint inputs, uint outputs ) {
        foreach( i ; 0..inputs ) {
            Node inputNode = addNode(Node.Type.input);
            foreach( o ; 0..outputs ) {
                Node outputNode = addNode(Node.Type.output);
                addConnection( i, o );
            }
        }
    }

    this( Genome g ) {
        foreach( n ; g.nodes ) {
            cloneNode(n);
        }
    }

    Genome dup() {
        return new Genome(this);
    }

    /// add node with type, without connections
    Node addNode( Node.Type type ) {
        Node n = new Node(type);
        nodes ~= n;
        return n;
    }

    /// add cloned node incl. connections
    Node cloneNode( Node n ) {
        Node newNode = new Node(n);
        nodes ~= newNode;

        return newNode;
    }

    Connection addConnection( uint startNode, uint endNode ) {
        Connection c = new Connection( startNode, endNode );
        nodes[startNode].addOutput(c);
        nodes[endNode].addInput(c);
        connections ~= c;
        return c;
    }

    Connection[] connections;
    Node[] nodes;
}