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
    float distance( Individuum i ) {
        float dx = (x - i.x);
        dx = dx*dx;
        float dy = (y - i.y);
        dy = dy*dy;
        return sqrt(dx+dy);
    }
}



void mainA() {
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

    //delegate bool cmp(i) { return distance(ind[0], i) < thresh;};
    while( !list.empty ) {
        writefln("Length: %s", list[].walkLength());
        auto range = list[];
        auto reference = range.front;
        while(!range.empty) {
            if( reference.distance(range.front) < thresh ) {
                if( sc >= species.length ) {
                    species.length++;
                }
                species[sc] ~= range.front;
                list.linearRemove(range.take(1));
            }
            range.popFront();
        }
        sc++;
    }
    foreach(aidx, a; species) {
        writeln(aidx);
        foreach(b; a) {
            writeln("  ", b.nr);
        }
    }


}