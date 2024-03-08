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