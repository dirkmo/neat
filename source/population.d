module neat.population;

import neat.genepool;
import neat.individual;

class Population {                                                           
    this( uint popsize, uint inputs, uint outputs, bool recurrent ) { 
        pool = new Genepool(
            inputs, outputs,
            true /*fullyConnected*/,
            recurrent /*enableRecurrentNets*/
        );  
        individuals.length = popsize;
        foreach( ref i; individuals ) { 
            i = new Individual( pool, true /*createConPhenotype*/ );
        }   
    }   
        
//private:

    Genepool pool;
    Individual[] individuals;
}