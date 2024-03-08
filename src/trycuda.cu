#include <stdio.h>
#include <assert.h>
#include <cuda.h>
#define Block_Size 256
#define NUM_BLOCKS_MAX 2147483647
__global__ void stretch(double *a, 
                        double *b,
                        int N){
    int index= blockIdx.x * blockDim.x + threadIdx.x;
    int i,j,stride2=N*N;
    double temp1,temp2;
    if(index<stride2){
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
    double adelt=0.0;
    int dij,drj,dis,drs;
    double temp,temp1;
    // if(index==0){printf("\n674374\n");}
    if(index<norb*stride3){
        int r = index / stride3;
        int i = (index / stride2) % norb;
        int s = (index / norb) % norb;
        int j = index % norb;
        dij = (i == j) ? 1.0 : 0.0;
        drj = (r == j) ? 1.0 : 0.0;
        dis = (i == s) ? 1.0 : 0.0;
        drs = (r == s) ? 1.0 : 0.0;
        adelt=0.0;
        adelt=aplus_ck[r*norb+s]*dij-aplus_ck[i*norb+s]*drj-aplus_ck[r*norb+j]*dis+
                    aplus_ck[i*norb+j]*drs;
        
        temp=hess_ck[r*stride3+i*stride2+s*norb+j];
        temp=temp+dh_ck[r*stride3+s*stride2+i*norb+j]-
        dh_ck[i*stride3+s*stride2+r*norb+j]-dh_ck[r*stride3+j*stride2+i*norb+s]+
        dh_ck[i*stride3+j*stride2+r*norb+s];
        temp=temp+
                Y_ck[r*stride3+i*stride2+s*norb+j]-Y_ck[i*stride3+r*stride2+s*norb+j]-
                Y_ck[r*stride3+i*stride2+j*norb+s]+Y_ck[i*stride3+r*stride2+j*norb+s];

        hess_ck[r*stride3+i*stride2+s*norb+j]=temp-0.5*adelt;

        if(s==norb&&j==norb){
            Hdiag_c[r*norb+i]=hess_ck[r*stride3+i*stride2+r*norb+i];
        }
    }
}
extern "C" void hess3c_(double *aplus_c,double *F_c,int *NP,
                        double *hess_c,double *dh_c,double *Y_c,
                        double *Hdiag_c){
    printf("\nEnter hess3c_ C++\n");
    int norb=*NP;
    int stride4=norb*norb*norb*norb,stride3=norb*norb*norb,stride2=norb*norb;
    int threads_per_block=Block_Size,num_blocks=(stride2+threads_per_block-1)/threads_per_block;
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

