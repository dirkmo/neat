module neat.connection;

import neat.node;

import std.random;
import std.string;

///
class ConGene {
    ///
    this( NodeGene start, NodeGene end ) {
        _innovation = _innovationCounter++;
        _startNode = start;
        _endNode = end;
    }

    ///
    NodeGene start()  { return _startNode; }
    ///
    NodeGene end()    { return _endNode; }

    uint innovation()        const @property { return _innovation; }

    override string toString() const {
        return format( "%s", _innovation );
    }

private:
    uint _innovation;
    NodeGene _startNode;
    NodeGene _endNode;
    
    static uint _innovationCounter;
}

///
class Connection {
    ///
    this( ConGene gene, Node input, Node output ) {
        this.gene = gene;
        _weight = uniform( -1.0f, 1.0f );
        enabled = true;
    }

    ///
    this( Connection con, Node input, Node output ) {
        this.gene = con.gene;
        _weight = con.weight();
        enabled = con.enabled;
        _input = input;
        _output = output;
    }

    const(Node) input()  const @property { return _input; }
    const(Node) output() const @property { return _output; }
    float weight()    const @property { return _weight; }
    uint innovation() const @property { return gene.innovation; }

    bool enabled;

private:
    ConGene gene;
    Node _input;
    Node _output;

    float _weight;
}