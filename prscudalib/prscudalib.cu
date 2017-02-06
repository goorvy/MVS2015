
#ifndef _PRS_CUDA_FUNC_
#define _PRS_CUDA_FUNC_

#include <dsp\bitbuf.h>
#include <vector>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#define imin(a,b) (a<b?a:b)

//__constant__ int seqs[16384]; // максимальный объем константной памяти на моей GTX750 = 64k

__global__ void akf_kernel(int *akf, const int *seq1, const int *seq2, const int seq_size) { // расчет сразу всей КФ (по всем сдвигам)
    //__shared__ int cache[THREADS_PER_BLOCK];
    extern __shared__ int cache[];
    int tid = threadIdx.x + blockIdx.y * blockDim.x; // 1-я размерность блоков расчитывает 1 значение КФ
    int cacheIndex = threadIdx.x;
    //int temp = seq1[tid] ^ seq2[(blockIdx.x + tid) % seq_size];
    int temp = seq1[tid] * seq2[(blockIdx.x + tid) % seq_size]; 
    cache[cacheIndex] = temp;
    __syncthreads();
    temp = blockDim.x / 2;
    while (temp != 0) {
        if (cacheIndex < temp)
            cache[cacheIndex] += cache[cacheIndex + temp];
        __syncthreads();
        temp /= 2;
    }
    if (cacheIndex == 0) // в 0вом элементе сумма всех элементов массива
        akf[blockIdx.x * gridDim.y + blockIdx.y] = cache[0]; // 2-я размерность блоков расчитывает сдвиги
}

namespace dsp {
    namespace prs {
        int cudaXcorr(std::vector<float> &akf,
            const dsp::BitBuffer<dsp::u32>  &seq1,
            const dsp::BitBuffer<dsp::u32>  &seq2,
            const bool                      normalize = true) {

            if (seq1.size() != seq2.size())
                return 128; // ошибка, последовательности разной длины

            int *dev_seq1 = 0;
            int *dev_seq2 = 0;
            int *dev_akf = 0;
            unsigned int threads_per_block = imin(1024, seq1.size());
            unsigned int blocks_per_seq = imin(65535, (seq1.size() + threads_per_block - 1) / threads_per_block);
            cudaError_t err;

            int *int_seq1 = new int[seq1.size()];
            int *int_seq2 = new int[seq1.size()];
            int *akf_temp = new int[seq1.size() * blocks_per_seq];

            for (int k = 0; k < seq1.size(); ++k) {
                //int_seq1[k] = seq1[k];
                //int_seq2[k] = seq2[k];
                if (seq1[k] == 1)
                    int_seq1[k] = 1;
                else
                    int_seq1[k] = -1;
                if (seq2[k] == 1)
                    int_seq2[k] = 1;
                else
                    int_seq2[k] = -1;
            }
            // GPU part.
            err = cudaMalloc((void**)&dev_seq1, seq1.size() * sizeof(int));
            //err = cudaMemcpyToSymbol(seqs, int_seq1, seq1.size() * sizeof(int));
            if (err != cudaSuccess) return (int)err;
            err = cudaMalloc((void**)&dev_seq2, seq1.size() * sizeof(int));
            if (err != cudaSuccess) return (int)err;
            err = cudaMalloc((void**)&dev_akf, seq1.size() * blocks_per_seq * sizeof(int));
            if (err != cudaSuccess) return (int)err;

            err = cudaMemcpy(dev_seq1, int_seq1, seq1.size() * sizeof(int), cudaMemcpyHostToDevice);
            if (err != cudaSuccess) return (int)err;
            err = cudaMemcpy(dev_seq2, int_seq2, seq1.size() * sizeof(int), cudaMemcpyHostToDevice);
            if (err != cudaSuccess) return (int)err;

            dim3 grids(seq1.size(), blocks_per_seq);
            dim3 threads(threads_per_block);
            akf_kernel <<< grids, threads, threads_per_block * sizeof(int) >>> (dev_akf, dev_seq1, dev_seq2, seq1.size());
            err = cudaGetLastError();
            if (err != cudaSuccess) return (int)err;

            err = cudaMemcpy(akf_temp, dev_akf, seq1.size() * blocks_per_seq * sizeof(int), cudaMemcpyDeviceToHost);
            if (err != cudaSuccess) return (int)err;

            for (int i = 0; i < seq1.size(); ++i) { // суммирование значений КФ по блокам в одно целое
                akf[i] = 0;
                for (int k = 0; k < blocks_per_seq; ++k)
                    akf[i] += akf_temp[i*blocks_per_seq + k];
            }

            if (normalize)
                for (int i = 0; i < seq1.size(); ++i)
                    akf[i] /= seq1.size();

            cudaFree(dev_akf);
            cudaFree(dev_seq1);
            cudaFree(dev_seq2);

            delete[] int_seq1;
            delete[] int_seq2;
            delete[] akf_temp;

            return 0; // нет ошибок
        };
    }
}

#endif

