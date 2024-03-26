#include <stdio.h>
#include <assert.h>
#include <cuda.h>
#define Block_Size 256
#define NUM_BLOCKS_MAX 2147483647

extern "C" void helloworldcuda_ (){
    printf("Hello world CUDA Routine !!\n");
    int deviceCount;
    int dev;
    cudaGetDeviceCount(&deviceCount);
    dev=0;
    cudaDeviceProp deviceProp;
    cudaGetDeviceProperties(&deviceProp, dev);
    printf("\nDevice %d: on %d \"%s\"\n", dev, deviceCount, deviceProp.name);
}
__global__ void kernel_aplus(double *a, double *b ,int N){
    int tid = blockIdx.x * blockDim.x + threadIdx.x;//j
    if(tid<N){
        if(tid==0){
            printf("67876\n");
        }
        for(int i=0;i<N;i++){
            a[tid*N+i]=b[tid*N+i]+b[i*N+tid];
        }
    }
    // for(int i=0;i<norb;i++){
    //     for(int j=0;j<norb;j++){
    //         // a[j][i]=F[j][i]+F[i][j]
    //         a[j*N+i]=b[j*N+i]+b[i*N+j]
    //     }
    // }
}
extern "C" void aplus_c_(double *a,double *b,int *NP){   
    printf("Enter GPU aplus  \n");
    int N=*NP;
    int size=N*N;
    int threads_per_block=Block_Size;
    int num_blocks=(N+threads_per_block-1)/threads_per_block;
    double *a_d,*b_d;
    cudaMalloc( (void **)&a_d, sizeof(double) * size );
    cudaMalloc( (void **)&b_d, sizeof(double) * size );

    cudaMemcpy( a_d, a, sizeof(double) * size, cudaMemcpyHostToDevice );
    cudaMemcpy( b_d, b, sizeof(double) * size, cudaMemcpyHostToDevice );
    
    kernel_aplus<<<num_blocks,threads_per_block>>>(a_d,b_d,N);
    cudaDeviceSynchronize();
    cudaMemcpy( a, a_d, sizeof(double) * size, cudaMemcpyDeviceToHost );
    
    cudaFree(a_d);
    cudaFree(b_d);
}

__device__ __host__ void mat_const12(double *GM1,int index_x,int index_y,double* T,double Dnum,int norb){
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
__global__ void kernel_dhgen(int *occ,int *tot,double *GM1,double *T,double *D,int nsub,int nact,int norb){
    int i = blockIdx.x * blockDim.x + threadIdx.x;//j
    int ioffset=0;
    int ioffset1=0;
    int i0,i1,i2,ix,joffset,joffset1,j0,j1,j2,jx;
    if(i<nsub){
        for(int ii=0;ii<occ[i];ii++){
            i0=ioffset;
            i1=tot[i];
            i2=occ[i];
            ix=ioffset1;
            joffset=0;
            joffset1=0;
            for(int j=0;j<nsub;j++){
                for(int jj=0;jj<occ[j];jj++){
                    j0=joffset;
                    j1=tot[j];
                    j2=occ[j];
                    jx=joffset1;
                    // GM1(i0+ii,j0+jj,:,:)=T*D(ix+ii,jx+jj)*2.0d0
                    mat_const12(GM1,(i0+ii),(j0+jj),T,D[(jx+jj)*nact+(ix+ii)],norb);
                }
            joffset=joffset+tot[i];
            joffset1=joffset1+occ[i];
            }
        }
        ioffset=ioffset+tot[i];
        ioffset1=ioffset1+occ[i];
    }
}
extern "C" void matgendh_(int *occ,int *tot,double *GM1,double *T,double *D,int *NP,int *NP2,int *NP3){
    printf("Enter GPU dhgen  \n");
    int nsub=*NP;
    int nact=*NP2;
    int norb=*NP3;
    int stride4=norb*norb*norb*norb;
    int stride2=norb*norb;
    int ioffset=0;
    int ioffset1=0;
    int i0,i1,i2,ix,joffset,joffset1,j0,j1,j2,jx=0;
    int threads_per_block=Block_Size;
    int num_blocks=(nsub+threads_per_block-1)/threads_per_block;
    int *d_occ,*d_tot;
    double *d_GM1,*d_T,*d_D;
    cudaMalloc( (void **)&d_occ, sizeof(int) * nsub );
    cudaMalloc( (void **)&d_tot, sizeof(int) * nsub );
    cudaMalloc( (void **)&d_GM1, sizeof(double) * stride4 );
    cudaMalloc( (void **)&d_T, sizeof(double) * stride2 );
    cudaMalloc( (void **)&d_D, sizeof(double) * nact*nact );
    cudaMemcpy( d_occ, occ, sizeof(int) * nsub, cudaMemcpyHostToDevice );
    cudaMemcpy( d_tot, tot, sizeof(int) * nsub, cudaMemcpyHostToDevice );
    cudaMemcpy( d_GM1, GM1, sizeof(double) * stride4, cudaMemcpyHostToDevice );
    cudaMemcpy( d_T, T, sizeof(double) * stride2, cudaMemcpyHostToDevice );
    cudaMemcpy( d_D, D, sizeof(double) * nact*nact, cudaMemcpyHostToDevice );
    kernel_dhgen<<<num_blocks,threads_per_block>>>(d_occ,d_tot,d_GM1,d_T,d_D,nsub,nact,norb);
    cudaDeviceSynchronize();
    cudaMemcpy( occ, d_occ, sizeof(int) * nsub, cudaMemcpyDeviceToHost );
    cudaMemcpy( tot, d_tot, sizeof(int) * nsub, cudaMemcpyDeviceToHost );
    cudaMemcpy( GM1, d_GM1, sizeof(double) * stride4, cudaMemcpyDeviceToHost );
    cudaMemcpy( T, d_T, sizeof(double) * stride2, cudaMemcpyDeviceToHost );
    cudaMemcpy( D, d_D, sizeof(double) * nact*nact, cudaMemcpyDeviceToHost );
    cudaFree(d_occ);
    cudaFree(d_tot);
    cudaFree(d_GM1);
    cudaFree(d_T);
    cudaFree(d_D);
}
extern "C" void matgendhcpu_(int *occ,int *tot,double *GM1,double *T,double *D,int *NP,int *NP2,int *NP3){
    int nsub=*NP;
    int nact=*NP2;
    int norb=*NP3;
    int ioffset=0;
    int ioffset1=0;
    int i0,i1,i2,ix,joffset,joffset1,j0,j1,j2,jx=0;
    for(int i=0;i<nsub;i++){
        for(int ii=0;ii<occ[i];ii++){
            i0=ioffset;
            i1=tot[i];
            i2=occ[i];
            ix=ioffset1;
            joffset=0;
            joffset1=0;
            for(int j=0;j<nsub;j++){
                for(int jj=0;jj<occ[j];jj++){
                    j0=joffset;
                    j1=tot[j];
                    j2=occ[j];
                    jx=joffset1;
                    // GM1(i0+ii,j0+jj,:,:)=T*D(ix+ii,jx+jj)*2.0d0
                    mat_const12(GM1,(i0+ii),(j0+jj),T,D[(jx+jj)*nact+(ix+ii)],norb);
                }
                joffset=joffset+tot[i];
                joffset1=joffset1+occ[i];
            }
        }
        ioffset=ioffset+tot[i];
        ioffset1=ioffset1+occ[i];
    }
}
void mat_TM1(double *TM1,double *U,double dtmp,int z,int q,int norb){
    int index_xy=z*norb*norb*norb+q*norb*norb;
    for(int i=0;i<norb;i++){
        for(int j=0;j<norb;j++){
            TM1[i*norb+j]+=U[index_xy+i*norb+j];
            // TM1=TM1+U(:,:,k0+kk,l0+ll)*dtmp
        }
    }
    
}
extern "C" void matgeny_(int *occ,int *tot,double *GM1,double *T,double *D,int *NP,int *NP2,
int *NP3,int *group,double *P,double *U){
    int nsub=*NP;
    int nact=*NP2;
    int norb=*NP3;
    int ioffset=0;
    int ioffset1=0;
    int i0,i1,i2,ix,joffset,joffset1,j0,j1,j2,jx,lx,l0,l1,l2,loffset,loffset1,kx,k0,k1,k2;
    int koffset,koffset1;   
    double dtmp;
    double *TM1,*TM2,*TM3;
    for(int i=0;i<nsub;i++){
        for(int ii=0;ii<occ[i];ii++){
            i0=ioffset;
            i1=tot[i];
            i2=occ[i];
            ix=ioffset1;
            joffset=0;
            joffset1=0;
            for(int j=0;j<nsub;j++){
                for(int jj=0;jj<occ[j];jj++){
                    j0=joffset;
                    j1=tot[j];
                    j2=occ[j];
                    jx=joffset1;

                    memset(TM1, 0, norb*norb * sizeof(double));
                    memset(TM2, 0, norb*norb * sizeof(double));
                    memset(TM3, 0, norb*norb * sizeof(double));

                    koffset=0;
                    koffset1=0;
                    for(int k=0;k<nsub;k++){
                        for(int kk=0;kk<occ[k];kk++){
                            k0=koffset;
                            k1=tot[k];
                            k2=occ[k];
                            kx=koffset1;

                            loffset=0;
                            loffset1=0;
                            for(int l=0;l<nsub;l++){
                                if(group[8*i+j]==group[8*k+l]){
                                    for(int ll=1;ll<occ[l];ll++){
                                        l0=loffset;
                                        l1=tot[l];
                                        l2=occ[l];
                                        lx=loffset1;
                                        dtmp=P[(ix+ii)*nact*nact*nact+(kx+kk)*nact*nact+(lx+ll)*nact+(jx+jj)];
                                        // TM1=TM1+U(:,:,k0+kk,l0+ll)*dtmp
                                    }
                                }
                                if(group[8*i+k]==group[8*j+l]){
                                    for(int ll=1;ll<occ[l];ll++){
                                        l0=loffset;
                                        l1=tot[l];
                                        l2=occ[l];
                                        lx=loffset1;
                                        dtmp=P[(ix+ii)*nact*nact*nact+(jx+jj)*nact*nact+(lx+ll)*nact+(kx+kk)]+
                                             P[(ix+ii)*nact*nact*nact+(lx+ll)*nact*nact+(jx+jj)*nact+(kx+kk)];
                                
                                        // TM2=TM2+U(:,k0+kk,l0+ll,:)*dtmp
                                    }
                                }
                                loffset=loffset+tot[l];
                                loffset1=loffset1+occ[l];
                            }
                        }
                        koffset=koffset+tot[k];
                        koffset1=koffset1+occ[k];
                    }
                    // GM1(i0+ii,:,j0+jj,:)=(TM1+TM2)*2.0d0
                    free(TM1);
                    free(TM2);
                    free(TM3);
                }
                joffset=joffset+tot[j];
                joffset1=joffset1+occ[j];
            }
        }
        ioffset=ioffset+tot[i];
        ioffset1=ioffset1+occ[i];
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
      ptrs[vector_lane][thread_lane] = Ap[row + thread_lane];
    }
    const int row_start = ptrs[vector_lane][0];  // same as row+start = Ap[row]
    const int row_end = ptrs[vector_lane][1];    // same as row_end = Ap[row+1]

    // initialize local sum
    double sum = 0.0;

    if (THREADS_PER_VECTOR == 32 && row_end - row_start > 32) {
      int jj = row_start - (row_start & (THREADS_PER_VECTOR - 1)) + thread_lane;
      if (jj >= row_start && jj < row_end) sum += Ax[jj] * x[Aj[jj]];
      for (jj += THREADS_PER_VECTOR; jj < row_end; jj += THREADS_PER_VECTOR) sum += Ax[jj] * x[Aj[jj]];
    } else {
      for (int jj = row_start + thread_lane; jj < row_end; jj += THREADS_PER_VECTOR) sum += Ax[jj] * x[Aj[jj]];
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
                                const double *Avals,
                                const int *Arows,
                                const int *Acols,
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
//   cudaMemcpy(y,b2,Anrows* sizeof(double),cudaMemcpyHostToDevice);//y

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


__global__ void kernel_dhgen(int *occ,int *tot,double *GM1,
                             double *T,double *D,int nsub,
                             int nact,int norb,int i){
    int ii = blockIdx.x * blockDim.x + threadIdx.x;//j
    int ioffset=0;
    int ioffset1=0;
    int i0,ix,joffset,joffset1,j0,jx;
    if(ii<occ[i]){
        i0=ioffset;
        ix=ioffset1;
        joffset=0;
        joffset1=0;
        for(int j=0;j<nsub;j++){
            for(int jj=0;jj<occ[j];jj++){
                j0=joffset;
                jx=joffset1;
                // GM1(i0+ii,j0+jj,:,:)=T*D(ix+ii,jx+jj)*2.0d0
                mat_const12(GM1,(i0+ii),(j0+jj),T,D[(jx+jj)*nact+(ix+ii)],norb);
            }
            joffset=joffset+tot[i];
            joffset1=joffset1+occ[i];
        }
    }
    ioffset=ioffset+tot[i];
    ioffset1=ioffset1+occ[i];
}

extern "C" void matgendh_(int *occ,int *tot,double *GM1,
                          double *T,double *D,int *NP,
                          int *NP2,int *NP3){
    printf("Enter GPU dhgen  \n");
    int nsub=*NP;
    int nact=*NP2;
    int norb=*NP3;
    int stride4=norb*norb*norb*norb;
    int stride2=norb*norb;
    int threads_per_block=THREADS_PER_BLOCK;
    int num_blocks=(nsub+threads_per_block-1)/threads_per_block;
    int *d_occ,*d_tot;
    double *d_GM1,*d_T,*d_D;
    cudaMalloc( (void **)&d_occ, sizeof(int) * nsub );
    cudaMalloc( (void **)&d_tot, sizeof(int) * nsub );
    cudaMalloc( (void **)&d_GM1, sizeof(double) * stride4 );
    cudaMalloc( (void **)&d_T, sizeof(double) * stride2 );
    cudaMalloc( (void **)&d_D, sizeof(double) * nact*nact );
    cudaMemcpy( d_occ, occ, sizeof(int) * nsub, cudaMemcpyHostToDevice );
    cudaMemcpy( d_tot, tot, sizeof(int) * nsub, cudaMemcpyHostToDevice );
    cudaMemcpy( d_GM1, GM1, sizeof(double) * stride4, cudaMemcpyHostToDevice );
    cudaMemcpy( d_T, T, sizeof(double) * stride2, cudaMemcpyHostToDevice );
    cudaMemcpy( d_D, D, sizeof(double) * nact*nact, cudaMemcpyHostToDevice );
    for (int  i = 0; i < nsub; i++)
    {
        kernel_dhgen<<<num_blocks,threads_per_block>>>(d_occ,d_tot,d_GM1,d_T,d_D,nsub,nact,norb,i);
    }
    cudaDeviceSynchronize();
    cudaMemcpy( occ, d_occ, sizeof(int) * nsub, cudaMemcpyDeviceToHost );
    cudaMemcpy( tot, d_tot, sizeof(int) * nsub, cudaMemcpyDeviceToHost );
    cudaMemcpy( GM1, d_GM1, sizeof(double) * stride4, cudaMemcpyDeviceToHost );
    cudaMemcpy( T, d_T, sizeof(double) * stride2, cudaMemcpyDeviceToHost );
    cudaMemcpy( D, d_D, sizeof(double) * nact*nact, cudaMemcpyDeviceToHost );
    cudaFree(d_occ);
    cudaFree(d_tot);
    cudaFree(d_GM1);
    cudaFree(d_T);
    cudaFree(d_D);
}

extern "C" void xjfserial_host_(int *rows, 
                            int *nz_num,
                            int * Ap, 
                            int * Aj, 
                            double *dddxjf, 
                            double *hfdjshx,    
                            double *y)    
{
    // printf("\nserial_host_\n");
    int num_rows=*rows;
    for (int i = 0; i < num_rows; i++){
        // printf("%lf ",hfdjshx[i]);
        const int row_start = Ap[i]-1;
        const int row_end   = Ap[i+1]-1;
        double sum = 0.0;
        for (int jj = row_start; jj < row_end; jj++) {            
            const int j = Aj[jj]-1;  //column index
            sum += hfdjshx[j] * dddxjf[jj];
        }
        y[i] = sum; 
    }
}
