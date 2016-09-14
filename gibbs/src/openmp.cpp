//
//  openmp.cpp
//  cpp_test
//
//  Created by yulele on 16/9/10.
//  Copyright (c) 2016年 yulele. All rights reserved.
//

#include <stdio.h>
#include <omp.h>

int main(int argc, char * argv[]) {
    #pragma omp parallel num_threads(5)
    {
        printf("Hello World!\n");
    }

    #pragma omp parallel for num_threads(4)
    for (int i = 0; i < 16; i ++)
        printf("thread_num=%d i=%d\n", omp_get_thread_num(), i);
    return 0;
}


