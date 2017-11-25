module neat.node;

import neat.connection;

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
    void addInputConGene( ConGene input ) { inputs ~= input; }

    ///
    void addOutputConGene( ConGene output ) { outputs ~= output; }

    ///
    ConGene[] getInputConGenes() { return inputs; }
    
    ///
    ConGene[] getOutputConGenes() { return outputs; }

    ///
    uint id() const @property { return _id; }
    
    ///
    int layerIndex() const @property { return _layerIndex; }
    
    ///
    int layerIndex( int index ) @property { return _layerIndex = index; }

    ///
    Type type() const @property { return _type; }

private:
    Type _type;
    uint _id;
    int _layerIndex;
    ConGene[] inputs;
    ConGene[] outputs;

    static uint _nodeId;
}

class Node {

    this( NodeGene gene ) {
        _gene = gene;
    }

    void addInputConnection( Connection input ) {
        inputs ~= input;
    }

    void addOutputConnection( Connection output ) {
        outputs ~= output;
    }

    const(NodeGene) gene() const @property { return _gene; }

private:
    NodeGene _gene;
    Connection[] inputs;
    Connection[] outputs;
}