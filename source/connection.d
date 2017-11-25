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

    uint innovation() const @property { return _innovation; }

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
        _gene = gene;
        _weight = uniform( -1.0f, 1.0f );
        enabled = true;
    }

    ///
    this( Connection con, Node input, Node output ) {
        _gene = con.gene;
        _weight = con.weight();
        enabled = con.enabled;
        _input = input;
        _output = output;
    }

    Node input()  @property { return _input; }
    Node output() @property { return _output; }
    ConGene gene() @property { return _gene; }
    float weight()    const @property { return _weight; }
    uint innovation() const @property { return _gene.innovation; }
    
    bool enabled;

private:
    ConGene _gene;
    Node _input;
    Node _output;

    float _weight;
}