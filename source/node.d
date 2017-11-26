module neat.node;

import neat.connection;

import std.algorithm;
import std.string;

///
class NodeGene {
    ///
    enum Type { input, hidden, output }

    ///
    this( Type type ) {
        this._id = _nodeId++;
        this._type = type;
        if( _type == Type.input ) {
            layerIndex = 0;
        } else {
            layerIndex = -1;
        }
    }

    ///
    void addInputConGene( ConGene input ) { _inputs ~= input; }

    ///
    void addOutputConGene( ConGene output ) { _outputs ~= output; }

    ///
    ConGene[] getInputConGenes() { return _inputs; }
    
    ///
    ConGene[] getOutputConGenes() { return _outputs; }

    ///
    bool isInputToNode( NodeGene ng, out ConGene con ) {
        auto icg = _outputs.find!( cg => cg.end().id == ng.id )();
        if( icg.length >0 ) {
            con = icg[0];
            return true;
        }
        return false;
    }

    ///
    uint id() const @property { return _id; }
    
    ///
    int layerIndex() const @property { return _layerIndex; }
    
    ///
    int layerIndex( int index ) @property { return _layerIndex = index; }

    ///
    Type type() const @property { return _type; }

    ///
    override string toString() const {
        return format("%s", _id);
    }

private:
    Type _type;
    uint _id;
    int _layerIndex;
    ConGene[] _inputs;
    ConGene[] _outputs;

    static uint _nodeId;
}

class Node {

    this( NodeGene gene ) {
        _gene = gene;
    }

    void addInputConnection( Connection input ) {
        _inputs ~= input;
    }

    void addOutputConnection( Connection output ) {
        _outputs ~= output;
    }

    ///
    bool isInputToNode( Node n, out Connection con ) {
        auto ic = _outputs.find!( c => c.end().id == n.id )();
        if( ic.length >0 ) {
            con = ic[0];
            return true;
        }
        return false;
    }

    NodeGene gene() @property { return _gene; }
    uint id() const @property { return _gene.id; }

    ///
    override string toString() const {
        return format("%s", _gene._id);
    }

    float value = 0.0f;

private:
    NodeGene _gene;
    Connection[] _inputs;
    Connection[] _outputs;
}