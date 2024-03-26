#include <stdio.h>
#include <assert.h>
#include <cuda.h>
#include <omp.h>
#include <sys/time.h>
#include <sys/stat.h>
#define THREADS_PER_BLOCK_AG 32
#define THREADS_PER_BLOCK_AG22 192
#define THREADS_PER_VECTOR 32

#define CUDA_CHECK_CALL(call)                                                                                     \
  {                                                                                                               \
    cudaError_t cudaStatus = call;                                                                                \
    if (cudaSuccess != cudaStatus)                                                                                \
      fprintf(stderr, "error CUDA RT call \"%s\" in line %d of file %s failed with  %s (%d).\n", #call, __LINE__, \
              __FILE__, cudaGetErrorString(cudaStatus), cudaStatus);                                              \
  }


__global__ void computeA_1(int *occ,int noffset,int noffset1,
                           int moffset,int moffset1,int joffset,int joffset1,
                           int n,int m,int j,int norb,int nocc,
                           double *A,double *T,double *D){

    const size_t VECTORS_PER_BLOCK = THREADS_PER_BLOCK_AG / THREADS_PER_VECTOR; //6 最好整除
    __shared__ volatile double sdata[VECTORS_PER_BLOCK * THREADS_PER_VECTOR + THREADS_PER_VECTOR / 2];

    const int thread_id = THREADS_PER_BLOCK_AG * blockIdx.x + threadIdx.x;  // global thread index
    const int thread_lane = threadIdx.x & (THREADS_PER_VECTOR - 1);      // thread index within the vector
    const int vector_id = thread_id / THREADS_PER_VECTOR;                // global vector index

    // const int vector_lane = threadIdx.x / THREADS_PER_VECTOR;            // vector index within the block
    // const int num_vectors = VECTORS_PER_BLOCK * gridDim.x;               // total number of active vectors

    int n1 = vector_id / occ[m];//row
    int m1 = vector_id % occ[m];//col

    double sum=0.0;
    if(THREADS_PER_VECTOR == 32 && occ[j] > 32){
        int jj = 0 - (0 & (THREADS_PER_VECTOR - 1)) + thread_lane;
        if (jj >= 0 && jj < occ[j]) sum += 2.0*T[(n1+noffset)+(jj+joffset)*norb]*D[(m1+moffset1)+(jj+joffset1)*nocc];
    }else{
        for(int jj=thread_lane;jj<occ[j];jj+=THREADS_PER_VECTOR){
            sum+=2.0*T[(n1+noffset)+(jj+joffset)*norb]*D[(m1+moffset1)+(jj+joffset1)*nocc];
        }
    }
    sdata[threadIdx.x] = sum;
    double tmp;
    if (THREADS_PER_VECTOR > 16) {tmp = sdata[threadIdx.x + 16];sum += tmp;sdata[threadIdx.x] = sum;}
    if (THREADS_PER_VECTOR > 8) {tmp = sdata[threadIdx.x + 8];sum += tmp;sdata[threadIdx.x] = sum;}
    if (THREADS_PER_VECTOR > 4) {tmp = sdata[threadIdx.x + 4];sum += tmp;sdata[threadIdx.x] = sum;}
    if (THREADS_PER_VECTOR > 2) {tmp = sdata[threadIdx.x + 2];sum += tmp;sdata[threadIdx.x] = sum;}
    if (THREADS_PER_VECTOR > 1) {tmp = sdata[threadIdx.x + 1];sum += tmp;sdata[threadIdx.x] = sum;}

    A[(n1+noffset)+(m1+moffset)*norb]+=sum;
}

__global__ void computeA_11(int *total,int *occ,int noffset,int noffset1,
                           int moffset,int moffset1,int joffset,int joffset1,
                           int n,int m,int j,int norb,int nocc,
                           double *A,double *T,double *D){
    int index= blockIdx.x * blockDim.x + threadIdx.x;
    int n1=index/occ[m];
    int m1=index%occ[m];
    if(index<occ[m]*total[n]){
        // printf("~~~%lf",D[index]);
        for(int j1=0;j1<occ[j];j1++){
            A[(n1+noffset)+(m1+moffset)*norb]=A[(n1+noffset)+(m1+moffset)*norb]
            +2.0*T[(n1+noffset)+(j1+joffset)*norb]*D[(m1+moffset1)+(j1+joffset1)*nocc];
        }
    }
    // for(int n1=0;n1<orb_total[n];n1++){
    //     for(int m1=0;m1<orb_occ[m];m1++){
    //         for(int j1=0;j1<orb_occ[j];j1++){
    //             A[(n1+noffset)+(m1+moffset)*norb]=A[(n1+noffset)+(m1+moffset)*norb]
    //             +2.0*T[(n1+noffset)+(j1+joffset)*norb]*D[(m1+moffset1)+(j1+joffset1)*nocc];
    //         }
    //     }
    // }
}

__global__ void computeA_22(int *total,int *occ,int noffset,int noffset1,
                           int moffset,int moffset1,int joffset,int joffset1,int koffset,int koffset1,int loffset,int loffset1,
                           int n,int m,int j,int k,int l,int norb,int nocc,
                           double *A,double *p,double *U,int stride2norb,int stride2nocc,int stride3nocc,int stride3norb){
    int index= blockIdx.x * blockDim.x + threadIdx.x;
    int n1=index/occ[m];
    int m1=index%occ[m];
    // if(index==10){
    
    // }
    if(index<occ[m]*total[n]){
        printf("~~~%lf",p[index]);
        for(int j1=0;j1<occ[j];j1++){
            for(int k1=0;k1<occ[k];k1++){
                for(int l1=0;l1<occ[l];l1++){
                    A[(n1+noffset)+(m1+moffset)*norb]=A[(n1+noffset)+(m1+moffset)*norb]
                    +2.0*p[(m1+moffset1)+(k1+koffset1)*nocc+(l1+loffset1)*stride2nocc+stride3nocc*(j1+joffset1)]
                    *U[(n1+noffset)+norb*(j1+joffset)+(k1+koffset)*stride2norb+(l1+loffset)*stride3norb];
                }
            }
        }
    }

    // for(int n1=0;n1<orb_total[n];n1++){
    //     for(int m1=0;m1<orb_occ[m];m1++){
    //         for(int j1=0;j1<orb_occ[j];j1++){
    //             for(int k1=0;k1<orb_occ[k];k1++){
    //                 for(int l1=0;l1<orb_occ[l];l1++){
    //                     A[(n1+noffset)+(m1+moffset)*norb]=A[(n1+noffset)+(m1+moffset)*norb]
    //                     +2.0*p[(m1+moffset1)+(k1+koffset1)*nocc+(l1+loffset1)*stride2nocc+stride3nocc*(j1+joffset1)]
    //                     *U[(n1+noffset)+norb*(j1+joffset)+(k1+koffset)*stride2norb+(l1+loffset)*stride3norb];
    //                 }
    //             }
    //         }
    //     }
    // }
    
}

extern "C" void operator_(int *nsub,double* A,double *T,double *D,
                          int *norbA,int *noccD,
                          int *orb_total,int *orb_occ,
                          double *p,double *U,double *G){
    timeval start,end,s1,e1;
    double runtime=0.0,r1=0.0;
    gettimeofday(&start,NULL);
    int orbnsub=*nsub;
    int norb=*norbA;
    int nocc=*noccD;
    int noffset=0,noffset1=0,moffset,moffset1,joffset,joffset1;
    int koffset=0,koffset1=0,loffset=0,loffset1=0;
    int stride3nocc=nocc*nocc*nocc,stride2nocc=nocc*nocc;
    int stride3norb=norb*norb*norb,stride2norb=norb*norb;
    int *occ_d,*total_d;
    double *A_d,*T_d,*D_d;
    CUDA_CHECK_CALL(cudaMalloc( (void **)&occ_d,sizeof(int)*orbnsub));
    CUDA_CHECK_CALL(cudaMalloc( (void **)&total_d,sizeof(int)*orbnsub));
    CUDA_CHECK_CALL(cudaMalloc( (void **)&A_d,sizeof(double)*stride2norb));
    CUDA_CHECK_CALL(cudaMalloc( (void **)&T_d,sizeof(double)*stride2norb));
    CUDA_CHECK_CALL(cudaMalloc( (void **)&D_d,sizeof(double)*stride2nocc));
    CUDA_CHECK_CALL(cudaMemcpy( occ_d,orb_occ,sizeof(int)*orbnsub,cudaMemcpyHostToDevice));
    CUDA_CHECK_CALL(cudaMemcpy( total_d,orb_total,sizeof(int)*orbnsub,cudaMemcpyHostToDevice));
    CUDA_CHECK_CALL(cudaMemcpy( A_d,A,sizeof(double)*stride2norb,cudaMemcpyHostToDevice));
    CUDA_CHECK_CALL(cudaMemcpy( T_d,T,sizeof(double)*stride2norb,cudaMemcpyHostToDevice));
    CUDA_CHECK_CALL(cudaMemcpy( D_d,D,sizeof(double)*stride2nocc,cudaMemcpyHostToDevice));
    for(int n=0;n<orbnsub;n++){
        moffset=0;
        moffset1=0;
        for(int m=0;m<orbnsub;m++){
            joffset=0;
            joffset1=0;
            for(int j=0;j<orbnsub;j++){
                // int block_num_AG=orb_total[n]*orb_occ[m];
                // computeA_1<<<block_num_AG,THREADS_PER_BLOCK_AG>>>(occ_d,noffset,noffset1,moffset,moffset1,joffset,joffset1,n,m,j,norb,nocc,A_d,T_d,D_d);
                int block_num_AG11=(orb_total[n]*orb_occ[m]+THREADS_PER_BLOCK_AG-1)/THREADS_PER_BLOCK_AG;
                computeA_11<<<block_num_AG11,THREADS_PER_BLOCK_AG>>>(total_d,occ_d,noffset,noffset1,moffset,moffset1,joffset,joffset1,n,m,j,norb,nocc,A_d,T_d,D_d);
                // for(int n1=0;n1<orb_total[n];n1++){
                //     for(int m1=0;m1<orb_occ[m];m1++){
                //         for(int j1=0;j1<orb_occ[j];j1++){
                //             A[(n1+noffset)+(m1+moffset)*norb]=A[(n1+noffset)+(m1+moffset)*norb]
                //             +2.0*T[(n1+noffset)+(j1+joffset)*norb]*D[(m1+moffset1)+(j1+joffset1)*nocc];
                //         }
                //     }
                // }
                cudaDeviceSynchronize();
                joffset=joffset+orb_total[j];
                joffset1=joffset1+orb_occ[j];
            }
            moffset=moffset+orb_total[m];
            moffset1=moffset1+orb_occ[m];
        }
        noffset=noffset+orb_total[n];
        noffset1=noffset1+orb_occ[n];
    }
    // cudaMemcpy( orb_occ,occ_d,sizeof(double)*orbnsub,cudaMemcpyDeviceToHost);
    cudaMemcpy( A,A_d,sizeof(double)*stride2norb,cudaMemcpyDeviceToHost);
    // cudaMemcpy( T,T_d,sizeof(double)*stride2norb,cudaMemcpyDeviceToHost);
    // cudaMemcpy( D,D_d,sizeof(double)*stride2nocc,cudaMemcpyDeviceToHost);
    cudaFree(occ_d);
    cudaFree(total_d);
    cudaFree(A_d);
    cudaFree(T_d);
    cudaFree(D_d);
    // double *p_d,*U_d;
    // CUDA_CHECK_CALL(cudaMalloc( (void **)&p_d,sizeof(double)*stride3nocc*nocc));
    // CUDA_CHECK_CALL(cudaMalloc( (void **)&U_d,sizeof(double)*stride3norb*norb));
    // CUDA_CHECK_CALL(cudaMemcpy(p_d,p,sizeof(double)*stride3nocc*nocc,cudaMemcpyHostToDevice));
    // CUDA_CHECK_CALL(cudaMemcpy(U_d,U,sizeof(double)*stride3norb*norb,cudaMemcpyHostToDevice));
    printf("***()*(*)(*&(*&*&^*&^^*&(*()*)&*(^*&^%&%&%$&$^%$%&(&))))");
    

    noffset=0;
    noffset1=0;
    moffset=0;moffset1=0;joffset=0;joffset1=0;
    for(int n=0;n<orbnsub;n++){
        moffset=0;
        moffset1=0;
        for(int m=0;m<orbnsub;m++){
            joffset=0;
            joffset1=0;
            for(int j=0;j<orbnsub;j++){
                koffset=0;
                koffset1=0;
                for(int k=0;k<orbnsub;k++){
                    loffset=0;
                    loffset1=0;
                    for(int l=0;l<orbnsub;l++){
                        // int block_num_AG22=(orb_total[n]*orb_occ[m]+THREADS_PER_BLOCK_AG-1)/THREADS_PER_BLOCK_AG;
                        // computeA_22<<<block_num_AG22,THREADS_PER_BLOCK_AG>>>(orb_total,orb_occ,noffset,noffset1,
                        //    moffset,moffset1,joffset,joffset1,koffset,koffset1,loffset,loffset1,
                        //    n,m,j,k,l,norb,nocc,
                        //    A_d,p_d,U_d,stride2norb,stride2nocc,stride3nocc,stride3norb);
                        // cudaDeviceSynchronize();
                        #pragma omp parallel for
                        for(int n1=0;n1<orb_total[n];n1++){
                            #pragma omp parallel for
                            for(int m1=0;m1<orb_occ[m];m1++){
                                for(int j1=0;j1<orb_occ[j];j1++){
                                    for(int k1=0;k1<orb_occ[k];k1++){
                                        for(int l1=0;l1<orb_occ[l];l1++){
                                            A[(n1+noffset)+(m1+moffset)*norb]=A[(n1+noffset)+(m1+moffset)*norb]
                                            +2.0*p[(m1+moffset1)+(k1+koffset1)*nocc+(l1+loffset1)*stride2nocc+stride3nocc*(j1+joffset1)]
                                            *U[(n1+noffset)+norb*(j1+joffset)+(k1+koffset)*stride2norb+(l1+loffset)*stride3norb];
                                        }
                                    }
                                }
                            }
                        }
                        loffset=loffset+orb_total[l];
                        loffset1=loffset1+orb_occ[l];
                    }
                    koffset=koffset+orb_total[k];
                    koffset1=koffset1+orb_occ[k];
                }
                joffset=joffset+orb_total[j];
                joffset1=joffset1+orb_occ[j];
            }
            moffset=moffset+orb_total[m];
            moffset1=moffset1+orb_occ[m];
        }
        noffset=noffset+orb_total[n];
        noffset1=noffset1+orb_occ[n];
    }
    // cudaMemcpy( A,A_d,sizeof(double)*stride2norb,cudaMemcpyDeviceToHost);
    cudaFree(occ_d);
    cudaFree(total_d);
    cudaFree(A_d);
    // cudaFree(p_d);
    // cudaFree(U_d);

    // for(int i=0;i<stride2norb;i++){
    //     printf("~~%lf",A[i]);
    // }

    int ioffset=0;
    for(int i=0;i<orbnsub;i++){
        #pragma omp parallel for
        for(int k=0;k<orb_total[i];k++){
            for(int l=0;l<orb_total[i];l++){
                if(fabs(A[(k + ioffset)+(l+ioffset)*norb]) < 1.0e-9){
                    A[(k + ioffset)+(l+ioffset)*norb]=0.0;
                }
            }
        }
        ioffset=ioffset+orb_total[i];
    }

    #pragma omp parallel for
    for(int i=0;i<norb;i++){
        for(int j=0;j<norb;j++){
            if(fabs(A[i+j*norb])<1.0e-9){
                A[i+j*norb]=0.0;
            }
        }
    }
    
    // gettimeofday(&start,NULL);
    noffset=0;
    noffset1=0;
    for(int n=0;n<orbnsub;n++){
        moffset=0;
        moffset1=0;
        for(int m=0;m<orbnsub;m++){
            joffset=0;
            joffset1=0;
            for(int j=0;j<orbnsub;j++){
                koffset=0;
                koffset1=0;
                for(int k=0;k<orbnsub;k++){
                    #pragma omp parallel for
                    for(int n1=0;n1<orb_total[n];n1++){
                        #pragma omp parallel for
                        for(int m1=0;m1<orb_occ[m];m1++){
                            #pragma omp parallel for
                            for(int j1=0;j1<orb_occ[j];j1++){
                                #pragma omp parallel for
                                for(int k1=0;k1<orb_total[k];k1++){
                                    int mr=m1+moffset1;
                                    int jr=j1+joffset1;
                                    int ni=n1+noffset;
                                    int mi=m1+moffset;
                                    int ji=j1+joffset;
                                    int ki=k1+koffset;
                                    G[ni+mi*norb+ji*stride2norb+ki*stride3norb]+=2.0*T[ni+ki*norb]*D[mr+jr*nocc];
                                }
                            }
                        }
                    }
                    koffset=koffset+orb_total[k];
                    koffset1=koffset1+orb_occ[k];
                }
                joffset=joffset+orb_total[j];
                joffset1=joffset1+orb_occ[j];
            }
            moffset=moffset+orb_total[m];
            moffset1=moffset1+orb_occ[m];
        }
        noffset=noffset+orb_total[n];
        noffset1=noffset1+orb_occ[n];
    }
    cudaDeviceSynchronize();
    gettimeofday(&end,NULL);
    runtime+=1e3*(end.tv_sec-start.tv_sec)+1e-3*(end.tv_usec-start.tv_usec);
    runtime/=1000;
    printf("G_1time:%lf s\n",runtime);


    gettimeofday(&s1,NULL);
    int kooffset=0,kooffset1=0;
    noffset=0;
    noffset1=0;
    for(int n=0;n<orbnsub;n++){
        moffset=0;
        moffset1=0;
        for(int m=0;m<orbnsub;m++){
            kooffset=0;
            kooffset1=0;
            for(int ko=0;ko<orbnsub;ko++){
                joffset=0;
                joffset1=0;
                for(int j=0;j<orbnsub;j++){
                    koffset=0;
                    koffset1=0;
                    for(int k=0;k<orbnsub;k++){
                        loffset=0;
                        loffset1=0;
                        for(int l=0;l<orbnsub;l++){
                            #pragma omp parallel for
                            for(int n1=0;n1<orb_total[n];n1++){
                                // #pragma omp parallel for
                                for(int m1=0;m1<orb_occ[m];m1++){
                                    for(int ko1=0;ko1<orb_occ[ko];ko1++){
                                        for(int j1=0;j1<orb_occ[j];j1++){
                                            for(int k1=0;k1<orb_occ[k];k1++){
                                                for(int l1=0;l1<orb_occ[l];l1++){
                                                    int mr=m1+moffset1;
                                                    int jr=j1+joffset1;
                                                    int kr=k1+koffset1;
                                                    int lr=l1+loffset1;
                                                    int ni=n1+noffset;
                                                    int mi=m1+moffset;
                                                    int koi=ko1+kooffset;
                                                    int ji=j1+joffset;
                                                    int ki=k1+koffset;
                                                    int li=l1+loffset;
                                                    G[ni+mi*norb+ji*stride2norb+koi*stride3norb]+=
                                                    2.0*p[mr+kr*nocc+lr*stride2nocc+jr*stride3nocc]
                                                    *U[ni+koi*norb+ki*stride2norb+li*stride3norb]
                                                    +2.0*2.0*p[mr+jr*nocc+lr*stride2nocc+kr*stride3nocc]
                                                    *U[ni+ki*norb+li*stride2norb+koi*stride3norb];
                                                    // printf("8 ");
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                            loffset=loffset+orb_total[l];
                            loffset1=loffset1+orb_occ[l];
                        }
                        koffset=koffset+orb_total[k];
                        koffset1=koffset1+orb_occ[k];
                    }
                joffset=joffset+orb_total[j];
                joffset1=joffset1+orb_occ[j];   
                }
                kooffset=kooffset+orb_total[ko];
                kooffset1=kooffset1+orb_occ[ko];
            }
            moffset=moffset+orb_total[m];
            moffset1=moffset1+orb_occ[m];
        }
        noffset=noffset+orb_total[n];
        noffset1=noffset1+orb_occ[n];
    }

    cudaDeviceSynchronize();
    gettimeofday(&e1,NULL);
    r1+=1e3*(e1.tv_sec-s1.tv_sec)+1e-3*(e1.tv_usec-s1.tv_usec);
    r1/=1000;
    printf("G_2time:%lf s\n",r1);
}



