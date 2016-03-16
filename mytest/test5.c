#include "phg.h"

#if (PHG_VERSION_MAJOR <= 0 && PHG_VERSION_MINOR < 9)
# undef ELEMENT
typedef SIMPLEX ELEMENT;
typedef USE_MPI 1;
#endif
 
#include <string.h>
#include <math.h>
int main(int argc, char *argv[])
{
    int m1, m2, m3, m4, m5;

    m1=-2;
    m2=-1;
    m3=0;
    m4=1;
    m5=2;

    m1=abs(m1);
    m2=abs(m2);
    m3=abs(m3);
    m4=abs(m4);
    m5=abs(m5);

    printf("m1=%d,m1%2=%d \n",m1,m1%2);
    printf("m2=%d,m2%2=%d \n",m2,m2%2);
    printf("m3=%d,m3%2=%d \n",m3,m3%2);
    printf("m4=%d,m4%2=%d \n",m4,m4%2);
    printf("m5=%d,m5%2=%d \n",m5,m5%2);

    printf("M_PI=%f \n",M_PI);

    int i=-1;
    i=-i;
    printf("i=%d \n",i);



    
    return 0;
}

