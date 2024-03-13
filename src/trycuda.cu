#include <stdio.h>
#include <assert.h>
#include <cuda.h>
#define THREADS_PER_BLOCK 256 
#define NUM_BLOCKS_MAX 2147483647 //2的31次方-1
#define THREADS_PER_VECTOR 32
#define HESSIAN3_SHARED_MEM_BLOCK_SIZE 6144


__global__ void stretch(double *a,
                        double *b,
                        int N){
    int index= blockIdx.x * blockDim.x + threadIdx.x;
    int totalThreads = blockDim.x * gridDim.x;
    int i,j,stride2=N*N;
    double temp1,temp2;
    for(;index<stride2;index+=totalThreads){
        i=index/N;
        j=index%N;
        temp1=b[index];
        temp2=b[j*N+i];
        a[index]=temp1+temp2;
    }
}

__global__ void kernel_update_hess_stretch(double *aplus_ck,int norb,
                                           double *hess_ck,double *dh_ck,double *Y_ck,
                                           double *Hdiag_c,int stride3,int stride2){
    int index= blockIdx.x * blockDim.x + threadIdx.x;
    int totalThreads = blockDim.x * gridDim.x;
    double adelt=0.0;
    int dij,drj,dis,drs;
    double temp,temp1,temp2;
    // if(index==0){printf("\n674374\n");}
    double dh_ck_rsij,dh_ck_isrj,dh_ck_rjis,dh_ck_ijrs;
    double Y_ck_risj,Y_ck_irsj,Y_ck_rijs,Y_ck_irjs;
    double aplus_ck_rs,aplus_ck_is,aplus_ck_rj,aplus_ck_ij;
    double hess_ck_riri;
    for(;index<norb*stride3;index+=totalThreads){
        int r = index / stride3;
        int i = (index / stride2) % norb;
        int s = (index / norb) % norb;
        int j = index % norb;
        dij = (i == j) ? 1.0 : 0.0;
        drj = (r == j) ? 1.0 : 0.0;
        dis = (i == s) ? 1.0 : 0.0;
        drs = (r == s) ? 1.0 : 0.0;
        adelt=0.0;
        aplus_ck_rs=aplus_ck[r*norb+s];
        aplus_ck_is=aplus_ck[i*norb+s];
        aplus_ck_rj=aplus_ck[r*norb+j];
        aplus_ck_ij=aplus_ck[i*norb+j];
        adelt=aplus_ck_rs*dij-aplus_ck_is*drj-aplus_ck_rj*dis+aplus_ck_ij*drs;
        temp=hess_ck[r*stride3+i*stride2+s*norb+j];
        dh_ck_rsij=dh_ck[r*stride3+s*stride2+i*norb+j];
        dh_ck_isrj=dh_ck[i*stride3+s*stride2+r*norb+j];
        dh_ck_rjis=dh_ck[r*stride3+j*stride2+i*norb+s];
        dh_ck_ijrs=dh_ck[i*stride3+j*stride2+r*norb+s];
        temp1=dh_ck_rsij-dh_ck_isrj-dh_ck_rjis+dh_ck_ijrs;
        Y_ck_risj=Y_ck[r*stride3+i*stride2+s*norb+j];
        Y_ck_irsj=Y_ck[i*stride3+r*stride2+s*norb+j];
        Y_ck_rijs=Y_ck[r*stride3+i*stride2+j*norb+s];
        Y_ck_irjs=Y_ck[i*stride3+r*stride2+j*norb+s];
        temp2=  Y_ck_risj-Y_ck_irsj-Y_ck_rijs+Y_ck_irjs;
        hess_ck[r*stride3+i*stride2+s*norb+j]=temp+temp1+temp2-0.5*adelt;
        if(s==norb&&j==norb){
            hess_ck_riri=hess_ck[r*stride3+i*stride2+r*norb+i];
            Hdiag_c[r*norb+i]=hess_ck_riri;
        }
    }
}

extern "C" void hess3c_(double *aplus_c,double *F_c,int *NP,
                        double *hess_c,double *dh_c,double *Y_c,
                        double *Hdiag_c){
    printf("\nEnter hess3c_ C++\n");
    int norb=*NP;
    int stride4=norb*norb*norb*norb,stride3=norb*norb*norb,stride2=norb*norb;
    int threads_per_block=THREADS_PER_BLOCK,num_blocks=(stride2+threads_per_block-1)/threads_per_block;
    num_blocks=min(NUM_BLOCKS_MAX,num_blocks);
    printf("num_blocks\n  %d ",num_blocks);
    double *aplus_c_d,*F_c_d;
    cudaMalloc( (void **)&aplus_c_d, sizeof(double) * stride2 );
    cudaMalloc( (void **)&F_c_d, sizeof(double) * stride2 );
    cudaMemcpy( aplus_c_d, aplus_c, sizeof(double) * stride2, cudaMemcpyHostToDevice );
    cudaMemcpy( F_c_d, F_c, sizeof(double) * stride2, cudaMemcpyHostToDevice );
    stretch<<<num_blocks,threads_per_block>>>(aplus_c_d,F_c_d,norb);
    cudaDeviceSynchronize();
    // for(int index=0;index<stride2;index++){
    //     int i=index/norb;
    //     int j=index%norb;
    //     aplus_c[index]=F_c[index]+F_c[j*norb+i];
    // }
    // for(int i=0;i<norb;i++){
    //     for(int j=0;j<norb;j++){
    //         aplus_c[j*norb+i]=F_c[j*norb+i]+F_c[i*norb+j];
    //     }
    // }
    
    double *hess_c_d,*dh_c_d,*Y_c_d,*Hdiag_c_d;
    cudaMalloc( (void **)&hess_c_d, sizeof(double) * stride4 );
    cudaMalloc( (void **)&dh_c_d, sizeof(double) * stride4 );
    cudaMalloc( (void **)&Y_c_d, sizeof(double) * stride4 );
    cudaMalloc( (void **)&Hdiag_c_d, sizeof(double) * stride2 );
    cudaMemcpy( hess_c_d, hess_c, sizeof(double) * stride4, cudaMemcpyHostToDevice );
    cudaMemcpy( dh_c_d, dh_c, sizeof(double) * stride4, cudaMemcpyHostToDevice );
    cudaMemcpy( Y_c_d, Y_c, sizeof(double) * stride4, cudaMemcpyHostToDevice );
    cudaMemcpy( Hdiag_c_d, Hdiag_c, sizeof(double) * stride2, cudaMemcpyHostToDevice );
    num_blocks=(stride4+threads_per_block-1)/threads_per_block;
    printf("num_block%d\n",num_blocks);
    kernel_update_hess_stretch<<<num_blocks,threads_per_block>>>(aplus_c_d,norb,
                                           hess_c_d,dh_c_d,Y_c_d,Hdiag_c_d,stride3,stride2);
    cudaDeviceSynchronize();
    
    cudaMemcpy( aplus_c, aplus_c_d, sizeof(double) * stride2, cudaMemcpyDeviceToHost );
    cudaMemcpy( hess_c, hess_c_d, sizeof(double) * stride4, cudaMemcpyDeviceToHost );
    cudaMemcpy( Hdiag_c, Hdiag_c_d, sizeof(double) * stride4, cudaMemcpyDeviceToHost );

    cudaFree(aplus_c_d);
    cudaFree(F_c_d);
    cudaFree(hess_c_d);
    cudaFree(dh_c_d);
    cudaFree(Y_c_d);
    cudaFree(Hdiag_c_d);

    // double adelt=0.0;
    // int dij=0,drj=0,dis=0,drs=0;
    // for(int index=0;index<norb*stride3;index++){
    //     int r = index / stride3;
    //     int i = (index / stride2) % norb;
    //     int s = (index / norb) % norb;
    //     int j = index % norb;
    //     dij = (i == j) ? 1.0 : 0.0;
    //     drj = (r == j) ? 1.0 : 0.0;
    //     dis = (i == s) ? 1.0 : 0.0;
    //     drs = (r == s) ? 1.0 : 0.0;
    //     adelt=0.0;
    //     adelt=aplus_c[r*norb+s]*dij-aplus_c[i*norb+s]*drj-aplus_c[r*norb+j]*dis+
    //                 aplus_c[i*norb+j]*drs;
        
    //     hess_c[r*stride3+i*stride2+s*norb+j]=
    //     hess_c[r*stride3+i*stride2+s*norb+j]+dh_c[r*stride3+s*stride2+i*norb+j]-
    //     dh_c[i*stride3+s*stride2+r*norb+j]-dh_c[r*stride3+j*stride2+i*norb+s]+
    //     dh_c[i*stride3+j*stride2+r*norb+s];

    //     hess_c[r*stride3+i*stride2+s*norb+j]=hess_c[r*stride3+i*stride2+s*norb+j]+
    //     Y_c[r*stride3+i*stride2+s*norb+j]-Y_c[i*stride3+r*stride2+s*norb+j]-
    //     Y_c[r*stride3+i*stride2+j*norb+s]+Y_c[i*stride3+r*stride2+j*norb+s];

    //     hess_c[r*stride3+i*stride2+s*norb+j]-=0.5*adelt;

    //     if(s==norb&&j==norb){
    //         Hdiag_c[r*norb+i]=hess_c[r*stride3+i*stride2+r*norb+i];
    //     }

    // }
    // for(int r=0;r<norb;r++){
    //     for(int i=0;i<norb;i++){
    //         for(int s=0;s<norb;s++){
    //             for(int j=0;j<norb;j++){
    //                 adelt=0.0;
    //                 dij = (i == j) ? 1.0 : 0.0;
    //                 drj = (r == j) ? 1.0 : 0.0;
    //                 dis = (i == s) ? 1.0 : 0.0;
    //                 drs = (r == s) ? 1.0 : 0.0;
    //                 adelt=aplus_c[r*norb+s]*dij-aplus_c[i*norb+s]*drj-aplus_c[r*norb+j]*dis+
    //                 aplus_c[i*norb+j]*drs;

    //                 hess_c[r*stride3+i*stride2+s*norb+j]=
    //                 hess_c[r*stride3+i*stride2+s*norb+j]+dh_c[r*stride3+s*stride2+i*norb+j]-
    //                 dh_c[i*stride3+s*stride2+r*norb+j]-dh_c[r*stride3+j*stride2+i*norb+s]+
    //                 dh_c[i*stride3+j*stride2+r*norb+s];

    //                 hess_c[r*stride3+i*stride2+s*norb+j]=hess_c[r*stride3+i*stride2+s*norb+j]+
    //                 Y_c[r*stride3+i*stride2+s*norb+j]-Y_c[i*stride3+r*stride2+s*norb+j]-
    //                 Y_c[r*stride3+i*stride2+j*norb+s]+Y_c[i*stride3+r*stride2+j*norb+s];

    //                 hess_c[r*stride3+i*stride2+s*norb+j]-=0.5*adelt;
    //             }
    //         }
    //         Hdiag_c[r*norb+i]=hess_c[r*stride3+i*stride2+r*norb+i];
    //     }
    // }
    printf("End C++\n");
}

__device__ __host__ void mat_const12(double *GM1,int index_x,
                                     int index_y,double* T,
                                     double Dnum,int norb){
    //T 是一个对称矩阵
    int index_xy=index_y*norb+index_x;
    int stride3=norb*norb*norb;
    int stride2=norb*norb;
    for(int i=0;i<norb;i++){
        for(int j=0;j<norb;j++){
            // int uu=j*stride3+i*stride2+index_xy;
            // printf("uu: %d,i: %d,j: %d\n",uu,i,j);
            GM1[j*stride3+i*stride2+index_xy]=T[j*norb+i]*Dnum*2.0;
            // GM1(i0+ii,j0+jj,:,:)=T*D(ix+ii,jx+jj)*2.0d0
            // printf("T[i*norb+j]:%lf  ,  T[i*norb+j]:%lf \n",T[i*norb+j],T[j*norb+i]);
        }
    }

}

__global__ void spmv_cusp(const unsigned int M,
                          const double *Ax,
                          const int *Ap,
                          const int *Aj,
                          const double *x,
                          double *y)
{
  const size_t VECTORS_PER_BLOCK = THREADS_PER_BLOCK / THREADS_PER_VECTOR;
  __shared__ volatile double sdata[VECTORS_PER_BLOCK * THREADS_PER_VECTOR + THREADS_PER_VECTOR / 2];
  __shared__ volatile int ptrs[VECTORS_PER_BLOCK][2];

  const int thread_id = THREADS_PER_BLOCK * blockIdx.x + threadIdx.x;  // global thread index
  const int thread_lane = threadIdx.x & (THREADS_PER_VECTOR - 1);      // thread index within the vector
  const int vector_id = thread_id / THREADS_PER_VECTOR;                // global vector index
  const int vector_lane = threadIdx.x / THREADS_PER_VECTOR;            // vector index within the block
  const int num_vectors = VECTORS_PER_BLOCK * gridDim.x;               // total number of active vectors

  for (int row = vector_id; row < M; row += num_vectors) {
    // use two threads to fetch Ap[row] and Ap[row+1],
    // considerably faster than the straightforward version
    if (thread_lane < 2) {
      ptrs[vector_lane][thread_lane] = Ap[row + thread_lane]-1;
    }
    const int row_start = ptrs[vector_lane][0];  // same as row+start = Ap[row]
    const int row_end = ptrs[vector_lane][1];    // same as row_end = Ap[row+1]

    // initialize local sum
    double sum = 0.0;

    if (THREADS_PER_VECTOR == 32 && row_end - row_start > 32) {
      int jj = row_start - (row_start & (THREADS_PER_VECTOR - 1)) + thread_lane;
      if (jj >= row_start && jj < row_end) sum += Ax[jj] * x[Aj[jj]-1];
      for (jj += THREADS_PER_VECTOR; jj < row_end; jj += THREADS_PER_VECTOR) sum += Ax[jj] * x[Aj[jj]-1];
    } else {
      for (int jj = row_start + thread_lane; jj < row_end; jj += THREADS_PER_VECTOR) sum += Ax[jj] * x[Aj[jj]-1];
    }
    // Store local sum in the shared memory
    sdata[threadIdx.x] = sum;
    // Reduce local sums to row sum
    double tmp;
    if (THREADS_PER_VECTOR > 16) {
      tmp = sdata[threadIdx.x + 16];
      sum += tmp;
      sdata[threadIdx.x] = sum;
    }
    if (THREADS_PER_VECTOR > 8) {
      tmp = sdata[threadIdx.x + 8];
      sum += tmp;
      sdata[threadIdx.x] = sum;
    }
    if (THREADS_PER_VECTOR > 4) {
      tmp = sdata[threadIdx.x + 4];
      sum += tmp;
      sdata[threadIdx.x] = sum;
    }
    if (THREADS_PER_VECTOR > 2) {
      tmp = sdata[threadIdx.x + 2];
      sum += tmp;
      sdata[threadIdx.x] = sum;
    }
    if (THREADS_PER_VECTOR > 1) {
      tmp = sdata[threadIdx.x + 1];
      sum += tmp;
      sdata[threadIdx.x] = sum;
    }
    // First thread writes the result
    if (thread_lane == 0) {y[row] = sdata[threadIdx.x];}
  }
}

extern "C"  void xjf_csr_(int *Anrowsc,
                                int *Annzc,
                                const int *Arows,
                                const int *Acols,
                                const double *Avals,
                                const double *xd,
                                double *b2)
{
  //定义A device
  int Anrows=*Anrowsc;
  int Annz=*Annzc;
  int *Arow_offset,*Aclo;
  double *Avalue,*xD,*y;

  cudaMalloc(((void **)(&xD)),Anrows* sizeof(double ));//
  cudaMemcpy(xD,xd,Anrows* sizeof(double ),cudaMemcpyHostToDevice);

  cudaMalloc(((void **)(&Avalue)),Annz* sizeof(double ));//
  cudaMemcpy(Avalue,Avals,Annz* sizeof(double ),cudaMemcpyHostToDevice);

  cudaMalloc(((void **)(&Arow_offset)),(Anrows+1)* sizeof(int ));//
  cudaMemcpy(Arow_offset,Arows,(Anrows+1)* sizeof(int ),cudaMemcpyHostToDevice);

  cudaMalloc(((void **)(&Aclo)),Annz* sizeof(int ));
  cudaMemcpy(Aclo,Acols,Annz* sizeof(int ),cudaMemcpyHostToDevice);//A.cols

  cudaMalloc(((void **)(&y)),Anrows* sizeof(double ));
  cudaMemcpy(y,b2,Anrows* sizeof(double),cudaMemcpyHostToDevice);//y

  const size_t VECTORS_PER_BLOCK  = THREADS_PER_BLOCK / THREADS_PER_VECTOR;//一个块中计算了多少行
  const size_t MAX_BLOCKS  = 2048;//cusp::system::cuda::detail::max_active_blocks
  const size_t NUM_BLOCKS = min(MAX_BLOCKS, (Anrows + (VECTORS_PER_BLOCK - 1)) / VECTORS_PER_BLOCK);
  
  spmv_cusp<<< NUM_BLOCKS,THREADS_PER_BLOCK,0 >>>(Anrows,Avalue,Arow_offset,Aclo,xD,y);
  cudaDeviceSynchronize();
  cudaMemcpy(b2,y,Anrows* sizeof(double ),cudaMemcpyDeviceToHost);
  
  cudaFree(Arow_offset);
  cudaFree(Aclo);
  cudaFree(Avalue);
  cudaFree(xD);
  cudaFree(y);
}

