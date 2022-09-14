#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <unistd.h>
#include <stdbool.h>
#include "bigbwt/gsa/gsacak.h"
#include "bigbwt/utils.h"


long n;         // length of Text[] (not including final 0 symbol)

Text = read_parse(arg.basename,&n);

// ------- compute largest input symbol (ie alphabet size-1)
long k=0;
for(long i=0;i<n;i++) {
if(Text[i]>k) k = Text[i];
}

// -------- alloc and compute SA of the parse
sa_index_t *SA = compute_SA(Text,n+1,k+1);

// transform SA->BWT inplace and write remapped last array, and possibly sainfo
sa_index_t *BWTsa = SA; // BWT overlapping SA
assert(n>1);
// first BWT symbol
assert(SA[0]==n);
BWTsa[0] = Text[n-1];

for(long i=1;i<=n;i++) {
    if(SA[i]==0) {
        assert(i==1);  // Text[0]=$abc... is the second lex word 
        BWTsa[i] = 0;   // eos in BWT, there is no phrase in D corresponding to this symbol so we write dummy values
    }
    else {
        BWTsa[i] = Text[SA[i]-1];
    }
}