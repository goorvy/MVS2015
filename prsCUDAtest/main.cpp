

#include <windows.h>
#include <prscudalib.h>
#include <conio.h>
#include <ctime>
#include <stdio.h>
#include <vector>
#include <dsp\bitbuf.h>
#include "C:\Users\dshubin\Documents\Visual Studio 2008\Projects\deBruijnFile\deBruijnFile\prs_debruijn_seqs.h"
#include <clocale>
#include <dsp\filebuf.h>

BOOL CALLBACK EnumWndProc(HWND hwnd, LPARAM lParam) { // ��� ��������� ����������� ���� �������
    if (GetWindowThreadProcessId(hwnd, NULL) == GetCurrentThreadId()) {
        *(HWND*)lParam = hwnd;
        return FALSE;
    }
    return TRUE;
}

#ifndef MAX_PATH
#define MAX_PATH 260
#endif

#define SEQ_LENGTH 2048
#define SEQ_INDEX 0

int main() {
    int seq_length;
    int seq_index;
    clock_t t_start, t_end;
    std::vector<float> akf;    

    std::setlocale(LC_CTYPE, "Russian_Russia.1251");

    printf("���������� ��� ����� ������������������ �� ������ �� GPU\n");

    printf("������� ����� ������������������ �� ������ = ");
    scanf("%d", &seq_length);

    printf("������� ������ ������������������ �� ������ = ");
    scanf("%d", &seq_index);

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
    dsp::BitBuffer<dsp::u32> bin_seq(seq_length);

    //printf("Sequence length = %d\nSequence index = %d\n", seq_length, seq_index);

    t_start = clock();
    dsp::prs::PRSDeBruijnSeq debruijn_seqs(seq_length);
    bin_seq = debruijn_seqs.get_seqs(seq_index);
    t_end = clock();

    printf("������������������ �� ������ ����� %d � �������� %d\n������������ �� %f ������.\n",
        seq_length, seq_index, ((float)t_end - (float)t_start) / CLOCKS_PER_SEC);

    t_start = clock();
    if (  int err = dsp::prs::cudaXcorr( akf, bin_seq, bin_seq, false )  ) {
        printf( "cudaError %d\n", err );
    }
    t_end = clock();

    float max = -seq_length;
    for (int i = 1; i < seq_length; ++i) {
        if (akf[i] > max)
            max = akf[i];
    }

    printf("Max AKF = %f\n", max / (float)seq_length);
    printf("������ ��� ���������� �� %f ������.\n", ((float)t_end - (float)t_start) / CLOCKS_PER_SEC);


    printf("\nFinish.");
    _getch();
    return 0;
}