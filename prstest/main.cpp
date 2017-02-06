
#include <conio.h>
#include <ctime>
#include <stdio.h>
#include <dsp\prs\PRSMath.h>
#include <vector>
#include <dsp\bitbuf.h>
#include "C:\Users\dshubin\Documents\Visual Studio 2008\Projects\deBruijnFile\deBruijnFile\prs_debruijn_seqs.h"
#include <clocale>


#define SEQ_LENGTH 2048
#define SEQ_INDEX 0

int main() {
    int seq_length;
    int seq_index;
    clock_t t_start, t_end;
    dsp::BitBuffer<dsp::u32> bin_seq;
    vector<int> akf;    

    std::setlocale(LC_CTYPE, "Russian_Russia.1251");

    printf("Вычисление АКФ одной последовательности де Брейна на CPU\n");

    printf("Введите длину последовательности де Брейна = ");
    scanf("%d", &seq_length);

    printf("Введите индекс последовательности де Брейна = ");
    scanf("%d", &seq_index);

    printf("Sequence length = %d\nSequence index = %d\n", seq_length, seq_index);

    akf.resize(seq_length);
    bin_seq.resize(seq_length);

    t_start = clock();
    dsp::prs::PRSDeBruijnSeq debruijn_seqs(seq_length);
    bin_seq = debruijn_seqs.get_seqs(seq_index);
    t_end = clock();

    printf("Последовательность де Брейна длины %d с индексом %d\nсформирована за %f секунд.\n",
        seq_length, seq_index, ((float)t_end - (float)t_start) / CLOCKS_PER_SEC);

    t_start = clock();
    akf = PRSMath::get_akf( &bin_seq );
    t_end = clock();

    float max = -seq_length;
    for (int i = 1; i < seq_length; ++i) {
        if (akf[i] > max)
            max = akf[i];
    }

    printf("Max AKF = %f\n", max / (float)seq_length);
    printf("Расчет АКФ произведен за %f секунд.\n", ((float)t_end - (float)t_start) / CLOCKS_PER_SEC);


    printf("\nFinish.");
    _getch();
    return 0;
}