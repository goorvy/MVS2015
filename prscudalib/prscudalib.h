#pragma once

#include <vector>
#include <dsp\bitbuf.h>

namespace dsp {
    namespace prs {
        int cudaXcorr(  std::vector<float>                &akf,             // �������� ������ �������� ���
                        const dsp::BitBuffer<dsp::u32>    &seq1,            // ������� ������������������ �1
                        const dsp::BitBuffer<dsp::u32>    &seq2,            // ������� ������������������ �2
                        const bool                        normalize = true);
    };
};
