#include <stdio.h>
#include <math.h>
#include <string.h>
#include "slave.h"

#include "div-arg.h"

#define J 64
#define I 1000

/*__thread_local volatile unsigned int get_reply,put_reply;*/

/*__thread_local int my_id;*/

__thread_local double a_slave[I],b_slave[I],c_slave[I];
__thread_local char buf_slave[1024];

extern double a[J][I],b[J][I],c[J][I];
extern unsigned long counter[64];

void mydiv(double *c, double *a, double *b) {
  *c = *a / *b;
}

void func(void *arg)
{
  int i,j;
  volatile unsigned int get_reply,put_reply;

  char *buf = (char *)arg;

  get_reply = 0;
  athread_get(PE_MODE,buf,buf_slave,sizeof(struct div1_t) + sizeof(struct div2_t),&get_reply,0,0,0);
  while(get_reply!=1);

  struct div1_t *pdvi1 = (struct div1_t *)buf;
  struct div2_t *pdiv2 = (struct div2_t *)(buf + sizeof(struct div1_t));
  struct div1_t *pdvi1_slave = (struct div1_t *)buf_slave;
  struct div2_t *pdiv2_slave = (struct div2_t *)(buf_slave + sizeof(struct div1_t));

  int my_id = athread_get_id(-1);
  if (my_id == 1) {
    printf("my_id: %d", my_id);
    printf("div1.foo1: %d, div1.bar1: %d\n", (int)pdvi1->foo1, pdvi1->bar1);
    printf("div2.foo2: %f, div2.bar2: %f\n", pdiv2->foo2, pdiv2->bar2);
    printf("div1.foo1_slave: %d, div1.bar1_slave: %d\n", (int)pdvi1_slave->foo1, pdvi1_slave->bar1);
    printf("div2.foo2_slave: %f, div2.bar2_slave: %f\n", pdiv2_slave->foo2, pdiv2_slave->bar2);
  }

  get_reply = 0;
  athread_get(PE_MODE,&a[my_id][0],&a_slave[0],I*8,&get_reply,0,0,0);
  athread_get(PE_MODE,&b[my_id][0],&b_slave[0],I*8,&get_reply,0,0,0);
  while(get_reply!=2);


  for(i=0;i<I;i++){
    /*c_slave[i]=a_slave[i]/b_slave[i];*/
    mydiv(&c_slave[i], &a_slave[i], &b_slave[i]);
  }

  put_reply=0;
  athread_put(PE_MODE,&c_slave[0],&c[my_id][0],I*8,&put_reply,0,0);
  while(put_reply!=1);
}
