#pragma once

#include <vector>
#include <dsp\bitbuf.h>

namespace dsp {
    namespace prs {
        int cudaXcorr(  std::vector<float>                &akf,             // выходной массив значений АКФ
                        const dsp::BitBuffer<dsp::u32>    &seq1,            // входная последовательность №1
                        const dsp::BitBuffer<dsp::u32>    &seq2,            // входная последовательность №2
                        const bool                        normalize = true);
    };
};
