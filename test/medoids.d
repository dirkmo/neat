module neat.medoids;

import std.algorithm;
import std.container;
import std.random;
import std.range;
import std.stdio;
import std.string;

///
interface Medoid {
    ///
    float distance( Medoid other );
    uint id() const @property;
}

///
class MedoidClassification {
    ///
    struct Cluster {
        ///
        Medoid medoid;
        ///
        Medoid[] members;
    }

    ///
    this( Medoid[] list, uint clusterCount, float thresh ) {
        this.list = list;
        this.clusterCount = clusterCount;
        this.thresh = thresh;
        // randomly pick medoids
        Medoid[] medoids;
        while( medoids.length < clusterCount ) {
            Medoid m = list[uniform(0, $)];
            if( !medoids.canFind!(l=>l.id == m.id)() ) {
                medoids ~= m;
            }
        }
        // create cluster structures
        foreach(m; medoids) {
            clusters ~= Cluster();
        }
        doClustering();
    }

    ///
    void doClustering() {
        // clear all clusters, keep medoids
        foreach(c; clusters) {
            c.members.length = 0;
        }
        // put items into cluster with nearest medoid
        foreach(it; list) {
            ulong minIdx;
            float minDist = float.max;
            foreach(idx, c; clusters) {
                const dist = it.distance(c.medoid);
                if( dist < minDist ) {
                    minIdx = idx;
                    minDist = dist;
                }
            }
            clusters[minIdx].members ~= it;
        }
        findNewMedoids();
    }

    ///
    void findNewMedoids() {
        foreach( cl; clusters ) {
            float minDist = float.max;
            Medoid newMediod;
            foreach( m; cl.members ) {
                const totalDist = calcTotalDistance( m, cl.members );
                if( totalDist < minDist ) {
                    newMediod = m;
                    minDist = totalDist;
                }
            }
            cl.medoid = newMediod;
        }
    }

    float calcTotalDistance( Medoid medoid, Medoid[] members ) {
        float dist = 0;
        foreach(m; members) {
            dist += medoid.distance(m);
        }
        return dist;
    }

    ///
    Cluster[] clusters;
    ///
    Medoid[] list;
    ///
    uint clusterCount;
    ///
    float thresh;
}

///
class Item : Medoid {
    ///
    this( float pos ) {
        this.pos = pos;
        _id = idCounter++;
    }

    ///
    float distance( Medoid m ) {
        Item it = cast(Item)m;
        return pos - it.pos;
    }

    ///
    uint id() const @property {
        return _id;
    }

private:
    float pos;
    uint _id;
    static uint idCounter = 0;
}

void main() {
    Medoid[] liste;
    liste.length = 20;
    write("Liste: ");
    foreach( ref l; liste ) {
        auto it = new Item(uniform(0, 100));
        l = it;
        writef("%s ", it.pos);
    }
    writeln();

    auto kmc = new MedoidClassification( liste, 5, 1.0f );
}

