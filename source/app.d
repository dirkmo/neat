import neat.genome;
import std.math;
import std.random;
import std.stdio;

float sigmoid( float x ) {
    return 1.0f / ( 1.0f + exp(-x) );
}

void printGenome(Genome g) {
	foreach( c; g.connections ) {
		writeln("net: ", c.innovation, ", i: ", c.inputNodeId, ", o: ", c.outputNodeId );
	}

	foreach( n; g.nodes ) {
		writeln("node: ", n.nodeId, ", type: ", n.type);
		write("  inputs: ");
		foreach( c; n.input ) {
			write(g.connections[c].innovation, " ");
		}
		write("\n  outputs: ");
		foreach( c; n.output ) {
			write(g.connections[c].innovation, " ");
		}
		writeln();
	}
}


void main()
{
	Genome g1 = new Genome(2, 1);
	writeln("g1:");

	printGenome(g1);
	
	Genome g2 = new Genome(g1);

	//g1.connections[0].innovation++;
	writeln("\ng2:");

	printGenome(g2);
	
	writeln("\ng3:");
	Genome g3 = g1.crossOver(g2, false);
	printGenome(g3);
}
