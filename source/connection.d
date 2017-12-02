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
        // add connection to nodes connection list
        _startNode.addOutputConGene(this);
        _endNode.addInputConGene(this);
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
    this( ConGene gene, Node start, Node end ) {
        _gene = gene;
        _start = start;
        _end = end;
        _weight = uniform( -1.0f, 1.0f );
        enabled = true;
        start.addOutputConnection(this);
        end.addInputConnection(this);
    }

    ///
    this( Connection con, Node start, Node end ) {
        _gene = con.gene;
        _weight = con.weight();
        _start = start;
        _end = end;
        enabled = con.enabled;
        start.addOutputConnection(this);
        end.addInputConnection(this);
    }

    float setWeight( float weight ) { return _weight = weight; }

    void mutateWeight( float probability, float strength ) {
        if( uniform(0.0f, 1.0f ) < probability ) {
            _weight += strength * uniform(-1.0f, 1.0f );
        }
    }

    override string toString() const {
        return format( "%s", innovation );
    }

    Node start()  @property { return _start; }
    Node end() @property { return _end; }
    ConGene gene() @property { return _gene; }
    float weight()    const @property { return _weight; }
    uint innovation() const @property { return _gene.innovation; }

    bool enabled;

private:
    ConGene _gene;
    Node _start;
    Node _end;

    float _weight;
}