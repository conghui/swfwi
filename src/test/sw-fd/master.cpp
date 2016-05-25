#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include "div-arg.h"

extern "C" {
#include "athread.h"

    void slave_func(void *);
}

static inline unsigned long rpcc()
{
    unsigned long time;
    asm("rtc %0": "=r" (time) : );
    return time;
}

#define J 64
#define I 1000


double a[J][I],b[J][I],c[J][I],cc[J][I];
double check[J];
unsigned long counter[J];

int main(void)
{
    int i,j;
    double checksum;
    double checksum2;
    unsigned long st,ed;

    printf("!!!!!!!!!! BEGIN INIT !!!!!!!!!!\n");fflush(NULL);
    for(j=0;j<J;j++)
        for(i=0;i<I;i++){
            a[j][i]=(i+j+0.5);
            b[j][i]=(i+j+1.0);
        }
    st=rpcc();
    for(j=0;j<J;j++)
        for(i=0;i<I;i++){
            cc[j][i]=(a[j][i])/(b[j][i]);
        }
    ed=rpcc();
    printf("the host     counter=%ld\n",ed-st);

    checksum=0.0;
    checksum2=0.0;
    athread_init();

    st=rpcc();

    struct div1_t div1;
    struct div2_t div2;
    div1.foo1 = 1;
    div1.bar1 = 2;
    div2.foo2 = 1.1;
    div2.bar2 = 2.2;

    char *buf = (char *)malloc(sizeof(div1) + sizeof(div2));
    memcpy(buf, &div1, sizeof(div1));
    memcpy(buf + sizeof(div1), &div2, sizeof(div2));

    //int np = 64;
    //for (int i = 0; i < np; i++) {
      //__real_athread_create(i, (void *)slave_func, 0);
    //}
    //for (int i = 0; i < np; i++) {
      //athread_wait(i);
    //}
    //for (int i = 0; i < np; i++) {
      //athread_end(i);
    //}
    __real_athread_spawn((void *)slave_func,buf);//fflush(NULL);
    athread_join();

    ed=rpcc();
    printf("the manycore counter=%ld\n",ed-st);

    printf("!!!!!!!!!! END JOIN !!!!!!!!!\n");fflush(NULL);

    for(j=0;j<J;j++)
        for(i=0;i<I;i++){
            checksum=checksum+c[j][i];
            checksum2=checksum2+cc[j][i];
        }

    printf("the master value is   %f!\n",checksum2);
    printf("the manycore value is %f!\n",checksum);

    if (fabs(checksum - checksum2) > 0.1) {
      printf("test error\n");
      exit(0);
    } else {
      printf("test pass!\n");
    }
    //athread_halt();
    printf("!!!!!!!!!! END HALT !!!!!!!!!\n");fflush(NULL);

}
