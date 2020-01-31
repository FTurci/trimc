#ifndef __SYS_H
#define __SYS_H

#include <omp.h>

template <class T> 
class A {
  public:

    A(){};
    int a;
    void compute();

};

template <class T> 
void A <T>::compute(){

  #pragma omp parallel num_threads(3)
  {
    const auto threadId = omp_get_thread_num();
    #pragma omp for 
    for(int i=0; i<10;i++) {
      printf("%d %d \n", omp_get_thread_num(), i);
    }
  }
}
#endif