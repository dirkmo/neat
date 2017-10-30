import neat.genome;
import std.math;
import std.random;
import std.stdio;

float sigmoid( float x ) {
    return 1.0f / ( 1.0f + exp(-x) );
}


void main()
{
	Genome g1 = new Genome(2, 2);
	writeln("g1:");
	foreach( c; g1.connections ) {
		writeln(c.innovation);
	}	
	
	Genome g2 = new Genome(g1);

	//g1.connections[0].innovation++;
	writeln("\ng2:");
	foreach( c; g2.connections ) {
		writeln(c.innovation);
	}
}
