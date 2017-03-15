#ifdef IN_R
#include "R.h"
#endif
#include "matrix.h"

int power (int base, int exp) {
    if (exp > 0) {
        return base * power(base, exp - 1);
    } else {
        return 1;
    }
}
