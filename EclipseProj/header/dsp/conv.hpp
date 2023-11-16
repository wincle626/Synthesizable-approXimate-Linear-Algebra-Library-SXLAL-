#ifndef CONV_CLION_CONV_H
#define CONV_CLION_CONV_H

template<class T, int L, int M>
void conv1d(T in[L], T out[L], T kernel[M])
{
    int i, j, k;

    // start convolution from out[M-1] to out[L-1] (last)
    for(i = M-1; i < L; ++i)
    {
        out[i] = 0;                             // init to 0 before accumulate

        for(j = i, k = 0; k < M; --j, ++k)
            out[i] += in[j] * kernel[k];
    }

    // convolution from out[0] to out[M-2]
    for(i = 0; i < M - 1; ++i)
    {
        out[i] = 0;                             // init to 0 before sum

        for(j = i, k = 0; j >= 0; --j, ++k)
            out[i] += in[j] * kernel[k];
    }

}

template<class T, int L, int M, int DX, int DY, int KX, int KY>
void conv2dslow(T in[L], T out[L], T kernel[M])
{
    int i, j, m, n, mm, nn;
    int kCenterX, kCenterY;                         // center index of kernel
    float sum;                                      // temp accumulation buffer
    int rowIndex, colIndex;

    // find center position of kernel (half of kernel size)
    kCenterX = KX / 2;
    kCenterY = KX / 2;

    for(i=0; i < DY; ++i)                // rows
    {
        for(j=0; j < DX; ++j)            // columns
        {
            sum = 0;                            // init to 0 before sum
            for(m=0; m < KY; ++m)      // kernel rows
            {
                mm = KY - 1 - m;       // row index of flipped kernel

                for(n=0; n < KX; ++n)  // kernel columns
                {
                    nn = KX - 1 - n;   // column index of flipped kernel

                    // index of input signal, used for checking boundary
                    rowIndex = i + (kCenterY - mm);
                    colIndex = j + (kCenterX - nn);

                    // ignore input samples which are out of bound
                    if(rowIndex >= 0 && rowIndex < DY && colIndex >= 0 && colIndex < DX)
                        sum += in[DX * rowIndex + colIndex] * kernel[KX * mm + nn];
                }
            }
            out[DX * i + j] = sum;
        }
    }

}

#endif //CONV_CLION_CONV_H
