module neat.medoids;

import std.algorithm;
import std.container;
import std.math;
import std.random;
import std.range;
import std.stdio;
import std.string;

///
interface Medoid(T) {
    float distance( Medoid m );
    uint id() const @property;    
}

/*
///
class Item : Medoid!Item {
    ///
    this( float pos ) {
        this.pos = pos;
        _id = idCounter++;
    }

    ///
    float distance( Medoid!Item m ) {
        Item it = cast(Item)m;
        return abs(pos - it.pos);
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
*/
///
class MedoidClassification(T) {
    ///
    class Cluster(T) {
        ///
        Medoid!T medoid;
        ///
        Medoid!T[] members;
        ///
        float meanDistance = float.max;
    }

    ///
    this( T[] list, uint clusterCount, float thresh ) {
        this.list = list;
        this.thresh = thresh;
        // randomly pick medoids
        T[] medoids;
        while( clusters.length < clusterCount ) {
            T m = list[uniform(0, $)];
            if( !medoids.canFind!(l=>l.id == m.id)() ) {
                auto newCluster = new Cluster!T();
                newCluster.medoid = m;
                clusters ~= newCluster;
                medoids ~= m;
            }
        }
        doClustering();
    }

    /// perform clustering until meanDistances is not improving anymode
    /// return: clusters
    T[][] clusterAll() {
        float dist = float.max;
        while( getMeanDistance() < dist ) {
            dist = getMeanDistance();
            doClustering();
        }
        T[][] clust;
        foreach(c; clusters) {
            T[] members;
            foreach( m; c.members) {
                members ~= cast(T)m;
            }
            clust ~= members;
        }
        return clust;
    }

    ///
    void doClustering() {
        // clear all clusters, keep medoids
        foreach(c; clusters) {
            c.members.length = 0;
        }
        // put Medoids into cluster with nearest medoid
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
            T newMediod;
            foreach( m; cl.members ) {
                const meanDist = calcMeanDistance( m, cl.members );
                if( meanDist < minDist ) {
                    newMediod = cast(T)m;
                    minDist = meanDist;
                }
            }
            cl.medoid = newMediod;
            cl.meanDistance = minDist;
        }
    }

    float calcMeanDistance( Medoid!T medoid, Medoid!T[] members ) {
        float dist = 0;
        foreach(m; members) {
            dist += medoid.distance(m);
        }
        return dist / members.length;
    }

    /// sum of mean distances of clusters
    float getMeanDistance() {
        float dist = 0;
        foreach(c; clusters) {
            dist += c.meanDistance;
        }
        return dist / clusters.length;
    }

    bool whichCluster( T item, out uint cluster ) {
        foreach( idx, c; clusters ) {
            if( c.members.canFind!(it=>it.id == item.id)() ) {
                cluster = cast(uint)idx;
                return true;
            }
        }
        return false;
    }

    T[] getCluster( uint idx ) {
        T[] items = cast(T[])clusters[idx].members;
        return items;
    }

    uint getClusterCount() @property { return cast(uint)clusters.length; }

    ///
    Cluster!T[] clusters;
    ///
    T[] list;
    ///
    float thresh;
}

/*
void main() {
    Item[] liste;
    liste.length = 100;
    write("Liste: ");
    foreach( ref l; liste ) {
        auto it = new Item(uniform(0.0f, 100.0f));
        l = it;
        writef("%s ", it.pos);
    }
    writeln();

    auto kmc = new MedoidClassification!Item( liste, 10, float.nan );
    float dist = float.max;
    uint i;
    while( kmc.getTotalDistance() < dist ) {
        dist = kmc.getTotalDistance();
        writeln("Durchgang ", i++);
        foreach(idx, c; kmc.clusters) {
            writef("Cluster %s, Medoid: %s: ", idx, (cast(Item)c.medoid).pos);
            foreach(m; c.members) {
                writef("%s ", (cast(Item)m).pos);
            }
            writeln();
        }
        kmc.doClustering();
        writeln("Totaldistance: ", kmc.getTotalDistance());
    }
}
*/
