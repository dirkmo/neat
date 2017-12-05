module neat.medoids;

import std.algorithm;
import std.container;
import std.math;
import std.random;
import std.range;
import std.stdio;
import std.string;

///
class Medoid {
    ///
    this( float pos ) {
        this.pos = pos;
        _id = idCounter++;
    }

    ///
    float distance( Medoid m ) {
        //Medoid it = cast(Medoid)m;
        return abs(pos - m.pos);
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

///
class MedoidClassification {
    ///
    class Cluster {
        ///
        Medoid medoid;
        ///
        Medoid[] members;
        ///
        float totalDistance = float.max;
    }

    ///
    this( Medoid[] list, uint clusterCount, float thresh ) {
        this.list = list;
        this.clusterCount = clusterCount;
        this.thresh = thresh;
        // randomly pick medoids
        Medoid[] medoids;
        while( clusters.length < clusterCount ) {
            Medoid m = list[uniform(0, $)];
            if( !medoids.canFind!(l=>l.id == m.id)() ) {
                auto newCluster = new Cluster();
                newCluster.medoid = m;
                clusters ~= newCluster;
                medoids ~= m;
            }
        }
        doClustering();
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
            Medoid newMediod;
            foreach( m; cl.members ) {
                const totalDist = calcTotalDistance( m, cl.members );
                if( totalDist < minDist ) {
                    newMediod = m;
                    minDist = totalDist;
                }
            }
            cl.medoid = newMediod;
            cl.totalDistance = minDist;
        }
    }

    float calcTotalDistance( Medoid medoid, Medoid[] members ) {
        float dist = 0;
        foreach(m; members) {
            dist += medoid.distance(m);
        }
        return dist;
    }

    float getTotalDistance() {
        float dist = 0;
        foreach(c; clusters) {
            dist += c.totalDistance;
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

void main() {
    Medoid[] liste;
    liste.length = 100;
    write("Liste: ");
    foreach( ref l; liste ) {
        auto it = new Medoid(uniform(0.0f, 100.0f));
        l = it;
        writef("%s ", it.pos);
    }
    writeln();

    auto kmc = new MedoidClassification( liste, 10, 1.0f );
    float dist = float.max;
    uint i;
    while( kmc.getTotalDistance() < dist ) {
        dist = kmc.getTotalDistance();
        writeln("Durchgang ", i++);
        foreach(idx, c; kmc.clusters) {
            writef("Cluster %s, Medoid: %s: ", idx, c.medoid.pos);
            foreach(m; c.members) {
                writef("%s ", m.pos);
            }
            writeln();
        }
        kmc.doClustering();
        writeln("Totaldistance: ", kmc.getTotalDistance());
    }
}

