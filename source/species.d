module neat.species;

import std.algorithm;
import std.container;
import std.math;
import std.random;
import std.range;
import std.stdio;

class Individuum {
    ulong nr;
    float x;
    float y;
}

float distance( Individuum i1, Individuum i2 ) {
    float dx = (i1.x - i2.x);
    dx = dx*dx;
    float dy = (i1.y - i2.y);
    dy = dy*dy;
    return sqrt(dx+dy);
}

void main() {
    Individuum[10] ind;
    foreach(idx, ref i; ind) {
        i = new Individuum();
        i.nr = idx;
        i.x = uniform(-10.0f, 10.0f);
        i.y = uniform(-10.0f, 10.0f);
        writeln(i);
    }
    float thresh = 10.0f;

    Individuum[][] species;
    uint sc;
    auto list = DList!Individuum(ind);

    while( true ) {

        auto range = list[].filter!( i => distance(ind[0], i) < thresh )();
        //writeln(range);
        if( sc >= species.length ) {
            species.length++;
        }
        foreach( r; range ) {
            species[sc] ~= r;
        }
        list.remove(range);
    }
}