#ifndef MERSENNE_H_
#define MERSENNE_H_

class MersenneTwister{
private:
  static const int N = 624;
  static const int M = 397;

  unsigned long mt[N]; /* the array for the state vector  */
  int mti; /* mti==N+1 means mt[N] is not initialized */

public:
  /* Initializing the array with current time.
     the seed is ignored.  */
  MersenneTwister();
  
  /* initialize by seed. */
  MersenneTwister(unsigned long seed);
  void initBySeed(unsigned long seed);

  /* Initialization by "sgenrand()" is an example. Theoretically,      */
  /* there are 2^19937-1 possible states as an intial state.           */
  /* This function allows to choose any of 2^19937-1 ones.             */
  /* Essential bits in "seed_array[]" is following 19937 bits:         */
  /*  (seed_array[0]&UPPER_MASK), seed_array[1], ..., seed_array[N-1]. */
  /* (seed_array[0]&LOWER_MASK) is discarded.                          */
  /* Theoretically,                                                    */
  /*  (seed_array[0]&UPPER_MASK), seed_array[1], ..., seed_array[N-1]  */
  /* can take any values except all zeros.                             */

  void lsgenrand(unsigned long seed_array[]);

  /* generating reals */
  double drand();
  double drand(double ma) {
    return ma * drand();
  }
  double drand(double lo, double hi) {
    return lo + (hi-lo) * drand();
  }

  /* generating ints */
  unsigned long irand();
  int irand(int ma);
  int irand(int lo, int hi);
};

#endif //MERSENNE_H_
