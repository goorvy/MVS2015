
//#include <windows.h>
//#include <dsp\filebuf.h>
#include <conio.h>
#include <ctime>
#include <stdio.h>
#include <vector>
#include <dsp\bitbuf.h>
#include "C:\Users\dshubin\Documents\Visual Studio 2008\Projects\deBruijnFile\deBruijnFile\prs_debruijn_seqs.h"
#include <clocale>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

//BOOL CALLBACK EnumWndProc(HWND hwnd, LPARAM lParam) { // ��� ��������� ����������� ���� �������
//    if (GetWindowThreadProcessId(hwnd, NULL) == GetCurrentThreadId()) {
//        *(HWND*)lParam = hwnd;
//        return FALSE;
//    }
//    return TRUE;
//}
//
//#ifndef MAX_PATH
//#define MAX_PATH 260
//#endif

#define imin(a,b) (a<b?a:b)

__global__ void akf_kernel(int *akf, const int *seq1, const int *seq2, const int seq_size) { // ������ ����� ���� �� (�� ���� �������)
                                                                                             //__shared__ int cache[THREADS_PER_BLOCK];
    extern __shared__ int cache[];
    int tid = threadIdx.x + blockIdx.y * blockDim.x; // 1-� ����������� ������ ����������� 1 �������� ��
    int cacheIndex = threadIdx.x;
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
    if (cacheIndex == 0) // � 0��� �������� ����� ���� ��������� �������
        akf[blockIdx.x * gridDim.y + blockIdx.y] = cache[0]; // 2-� ����������� ������ ����������� ������
}

int main()
{
    int seq_length;
    int seq_index;
    int seqs_amount;
    clock_t t_start, t_end;
    std::vector<float> akf;

    std::setlocale(LC_CTYPE, "Russian_Russia.1251");

    printf("���������� ��� ����� ������������������ �� ������ �� GPU\n");

    printf("������� ����� ������������������ �� ������ = ");
    scanf("%d", &seq_length);

    printf("������� ��������� ������ ������������������� �� ������ = ");
    scanf("%d", &seq_index);

    printf("������� ���������� ������������������� �� ������ = ");
    scanf("%d", &seqs_amount);

    //OPENFILENAME OpenFileName; // ��������� ��� �������
    //static TCHAR openfilename[255]; // ����� ����� �����
    //HWND hWnd; // ���������� ����
    //           // ����� ����� ����� ��� ����������
    //EnumWindows(EnumWndProc, (LPARAM)&hWnd); // ��������� ����������� ����������� ����
    //ZeroMemory(&OpenFileName, sizeof(OPENFILENAME));
    //OpenFileName.lStructSize = OPENFILENAME_SIZE_VERSION_400A;
    //OpenFileName.hwndOwner = hWnd;
    //OpenFileName.lpstrFile = openfilename;
    //OpenFileName.nMaxFile = MAX_PATH;
    //OpenFileName.lpstrFilter = "Binary Files\0*.bin\0\0"; // ������ ���� ������ � �������
    //OpenFileName.nFilterIndex = 1;
    //OpenFileName.lpstrFileTitle = NULL;
    //OpenFileName.nMaxFileTitle = 0;
    //OpenFileName.lpstrInitialDir = "C:\\Users\\dshubin\\Documents\\�����\\��� �� ������\\Save";
    //OpenFileName.Flags = OFN_PATHMUSTEXIST | OFN_FILEMUSTEXIST | OFN_NOCHANGEDIR;

    //if (!GetOpenFileName(&OpenFileName)) { // ����� ������� ������ ����� ��� ���������� �����
    //    cout << "������ �� ����� ��������.\n"; // � ������� ������ "������" �������������� ������������ ��������� �� ������ � ����
    //    return 1;
    //}
    //dsp::u32IFBB ifb(OpenFileName.lpstrFile);
    //u32* buf = new u32[SEQ_LENGTH];
    //for (int i = 0; i < SEQ_LENGTH; ++i)
    //    buf[i] = 0;
    //ifb(buf, SEQ_LENGTH);
    //dsp::BitBuffer<dsp::u32> bin_seq(buf, SEQ_LENGTH);

    akf.resize(seq_length);
    dsp::BitBuffer<dsp::u32> bin_seq1(seq_length);
    dsp::BitBuffer<dsp::u32> bin_seq2(seq_length);

    //printf("Sequence length = %d\nSequence index = %d\n", seq_length, seq_index);

    t_start = clock();
    dsp::prs::PRSDeBruijnSeq debruijn_seqs(seq_length);
    bin_seq1 = debruijn_seqs.get_seqs(seq_index);
    t_end = clock();

    printf("������������������ �� ������ ����� %d � �������� %d\n������������ �� %f ������.\n",
        seq_length, seq_index, ((float)t_end - (float)t_start) / CLOCKS_PER_SEC);


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
        akf_kernel << < grids, threads, threads_per_block * sizeof(int) >> > (dev_akf, dev_seq1, dev_seq2, seq1.size());
        err = cudaGetLastError();
        if (err != cudaSuccess) return (int)err;

        err = cudaMemcpy(akf_temp, dev_akf, seq1.size() * blocks_per_seq * sizeof(int), cudaMemcpyDeviceToHost);
        if (err != cudaSuccess) return (int)err;

        for (int i = 0; i < seq1.size(); ++i) { // ������������ �������� �� �� ������ � ���� �����
            akf[i] = 0;
            for (int k = 0; k < blocks_per_seq; ++k)
                akf[i] += akf_temp[i*blocks_per_seq + k];
        }

        cudaFree(dev_akf);
        cudaFree(dev_seq1);
        cudaFree(dev_seq2);

        delete[] int_seq1;
        delete[] int_seq2;
        delete[] akf_temp;

    return 0;
}

