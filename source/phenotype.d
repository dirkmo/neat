module neat.phenotype;

import neat.connection;
import neat.genepool;
import neat.node;

import std.algorithm;
import std.container;
import std.random;
import std.range;
import std.stdio;

class Phenotype {

    this() @disable {};

    this( Genepool pool, bool createConnections ) {
        this.pool = pool;
        if( createConnections ) {
            foreach( cg; pool.getConGenes() ) {
                addConnectionFromGenes(cg, cg.start, cg.end);
            }
        }
        // sort nodes by id, making input nodes first, followed by
        // output nodes
        nodes.sort!"a.id<b.id"();
    }

    /// create an identical copy of this
    Phenotype clone() {
        Phenotype newp = new Phenotype(pool, false);
        foreach( c; cons ) {
            auto nid1 = c.start.id;
            auto nid2 = c.end.id;

            Node n1new, n2new;
            auto narr = nodes.find!(n=>n.id == nid1)();
            if ( narr.length == 0 ) {
                n1new = new Node(c.start.gene);
            } else {
                n1new = narr[0];
            }
            narr = nodes.find!(n=>n.id == nid2)();
            if( narr.length == 0 ) {
                n2new = new Node(c.end.gene);
            }
            auto cg = new Connection( c.gene, n1new, n2new );
            newp.nodes ~= [n1new, n2new];
            newp.cons ~= cg;
        }

        return newp;
    }

    /// perform split up mutation
    void mutateSplitUpConnection(Connection con) {
        writeln(__FUNCTION__);
        if( !con.enabled ) {
            // gene is disabled, is already split up
            return;
        }
        ConGene cg1, cg2;
        pool.mutateSplitUpConGene(con.gene, cg1, cg2);
        // disable old connection
        con.enabled = false;
        auto ng = cg1.end;
        assert( !nodes.canFind!(n=>(n.id==ng.id))() );
        // add new node
        auto n3 = new Node(ng);
        nodes ~= n3;
        // add new connections
        auto n1 = con.start;
        auto n2 = con.end;
        auto c1 = new Connection( cg1, n1, n3 );
        auto c2 = new Connection( cg2, n3, n2 );
        cons ~= [c1, c2];
        c1.setWeight(1);
        c2.setWeight(con.weight);
    }

    ///
    void mutateAddConnection(Node n1, Node n2) {
        writeln(__FUNCTION__);
        Connection con;
        if( n1.isInputToNode( n2, con ) ) {
            // nodes already connected
            writefln("nodes %s and %s already connected.", n1, n2);
            return;
        }
        ConGene newConGene;
        if( pool.mutateAddNewConGene(n1.gene, n2.gene, newConGene) ) {
            // create connection from gene
            con = new Connection(newConGene, n1, n2);
            cons ~= con;
        }
    }

    /// mutate weights by probability with strength
    void mutateWeight(float probability, float strength) {
        // go over all connection phenotypes
        foreach(c; cons) {
            c.mutateWeight(probability, strength);
        }
    }

    ///
    Phenotype crossOver( Phenotype p2 ) {
        writeln(__FUNCTION__);
        Phenotype offspring = new Phenotype(pool, false);
        // add connections of p2 to slist
        auto lcp2 = SList!Connection(p2.cons);
        // loop over all connections of parent p1 (this) and look,
        // if p2 has same the connection gene
        foreach( i; cons ) {
            auto range = lcp2[].find!(c => c.innovation==i.innovation)();
            Connection c;
            if( range.empty ) {
                writeln("Only p1: ", i.innovation);
                // only p1 has the gene, take it
                c = i;
            } else {
                // both parents p1 and p2 have the gene
                // pick one at random
                writefln("both: %s", i.innovation);
                c = uniform(0.0f, 1.0f) < 0.5f ? i : range.front;
                // and remove connection from list lcp2
                lcp2.linearRemove(range.take(1));
                writeln("Length: ", lcp2.array().length);
            }
            // add connection (and possibly nodes) to offspring
            offspring.addConnection(i, i.start.gene, i.end.gene);
        }
        
        // the remaining connections in lcp2 are only in p2.
        // offspring will inherit them all
        foreach( i; lcp2 ) {
            // add connection (and possibly nodes) to offspring
            offspring.addConnection(i, i.start.gene, i.end.gene);
            writeln("Only p2: ", i.innovation);
        }
        writeln("cons.count = ", cons.length);
        return offspring;
    }

    private Connection addConnectionFromGenes(ConGene cg, NodeGene ng1, NodeGene ng2) {
        Node n1, n2;
        // does a node with ng1 gene exist?
        auto range = nodes.find!(n=>n.id==ng1.id)();
        if( range.empty ) {
            // no node with ng1 gene exists. create and add it.
            n1 = new Node(cg.start);
            nodes ~= n1;
        } else {
            // node with ng1 gene exists.
            n1 = range.front;
        }
        // does a node with ng2 gene exist?
        range = nodes.find!(n=>n.id==ng2.id)();
        if( range.empty ) {
            // no node with ng2 gene exists. create and add it.
            n2 = new Node(cg.end);
            nodes ~= n2;
        } else {
            // node with ng2 gene exists.
            n2 = range.front;
        }
        // finally add connection
        Connection newcon = new Connection(cg, n1, n2);
        cons ~= newcon;
        return newcon;
    }

    /// Add a copy of connection c, create nodes when necessary
    private void addConnection( Connection c, NodeGene ng1, NodeGene ng2 ) {
        writeln(__FUNCTION__);
        Connection con = addConnectionFromGenes( c.gene, ng1, ng2);
        con.setWeight(c.weight);
        con.enabled = c.enabled;
    }

//protected:
    Genepool pool;
    Connection[] cons;
    Node[] nodes;
}