#include <algorithm>    // std::fill_n

namespace mikan{

template<typename T>
T* create_1d_array(int len1)
{
    T *a1 = new T [len1];
    std::fill_n(a1, len1, (T)0);
    return a1;
}

template<typename T>
T** create_2d_array(int len2, int len1)
{
    T **a2 = new T* [len2];
    for (int i = 0; i < len2; ++i)
    {
        a2[i] = create_1d_array<T>(len1);
    }

    return a2;
}

template<typename T>
T*** create_3d_array(int len3, int len2, int len1)
{
    T ***a3 = new T ** [len3];
    for (int i = 0; i < len3; ++i)
    {
        a3[i] = create_2d_array<T>(len2, len1);
    }

    return a3;
}

template<typename T>
T**** create_4d_array(int len4, int len3, int len2, int len1)
{
    T ****a4 = new T *** [len4];
    for (int i = 0; i < len4; ++i)
    {
        a4[i] = create_3d_array<T>(len3, len2, len1);
    }

    return a4;
}

template<typename T>
T***** create_5d_array(int len5, int len4, int len3, int len2, int len1)
{
    T *****a5 = new T **** [len5];
    for (int i = 0; i < len5; ++i)
    {
        a5[i] = create_4d_array<T>(len4, len3, len2, len1);
    }

    return a5;
}

template<typename T>
T****** create_6d_array(int len6, int len5, int len4, int len3, int len2, int len1)
{
    T ******a6 = new T ***** [len6];
    for (int i = 0; i < len6; ++i)
    {
        a6[i] = create_5d_array<T>(len5, len4, len3, len2, len1);
    }

    return a6;
}

template<typename T>
T******* create_7d_array(int len7, int len6, int len5, int len4, int len3, int len2, int len1)
{
    T *******a7 = new T ****** [len7];
    for (int i = 0; i < len7; ++i)
    {
        a7[i] = create_6d_array<T>(len6, len5, len4, len3, len2, len1);
    }

    return a7;
}

template<typename T>
T******** create_8d_array(int len8, int len7, int len6, int len5, int len4, int len3, int len2, int len1)
{
    T ********a8 = new T ******* [len8];
    for (int i = 0; i < len8; ++i)
    {
        a8[i] = create_7d_array<T>(len7, len6, len5, len4, len3, len2, len1);
    }

    return a8;
}

template<typename T>
void delete_1d_array(T *a1)
{
    delete [] a1;
}

template<typename T>
void delete_2d_array(T **a2, int len2)
{
    for (int i = 0; i < len2; ++i)
    {
        delete_1d_array<T>(a2[i]);
    }

    delete [] a2;
}

template<typename T>
void delete_3d_array(T ***a3, int len3, int len2)
{
    for (int i = 0; i < len3; ++i)
    {
        delete_2d_array<T>(a3[i], len2);
    }

    delete [] a3;
}

template<typename T>
void delete_4d_array(T ****a4, int len4, int len3, int len2)
{
    for (int i = 0; i < len4; ++i)
    {
        delete_3d_array<T>(a4[i], len3, len2);
    }

    delete [] a4;
}

template<typename T>
void delete_5d_array(T *****a5, int len5, int len4, int len3, int len2)
{
    for (int i = 0; i < len5; ++i)
    {
        delete_4d_array<T>(a5[i], len4, len3, len2);
    }

    delete [] a5;
}

template<typename T>
void delete_6d_array(T ******a6, int len6, int len5, int len4, int len3, int len2)
{
    for (int i = 0; i < len6; ++i)
    {
        delete_5d_array<T>(a6[i], len5, len4, len3, len2);
    }

    delete [] a6;
}

template<typename T>
void delete_7d_array(T *******a7, int len7, int len6, int len5, int len4, int len3, int len2)
{
    for (int i = 0; i < len7; ++i)
    {
        delete_6d_array<T>(a7[i], len6, len5, len4, len3, len2);
    }

    delete [] a7;
}

template<typename T>
void delete_8d_array(T ********a8, int len8, int len7, int len6, int len5, int len4, int len3, int len2)
{
    for (int i = 0; i < len8; ++i)
    {
        delete_7d_array<T>(a8[i], len7, len6, len5, len4, len3, len2);
    }

    delete [] a8;
}

// Explicit template instantiation
template float* create_1d_array<float>(int len1);
template float** create_2d_array<float>(int len2, int len1);
template float*** create_3d_array<float>(int len3, int len2, int len1);
template float**** create_4d_array<float>(int len4, int len3, int len2, int len1);
template float***** create_5d_array<float>(int len5, int len4, int len3, int len2, int len1);
template float****** create_6d_array<float>(int len6, int len5, int len4, int len3, int len2, int len1);
template float******* create_7d_array<float>(int len7, int len6, int len5, int len4, int len3, int len2, int len1);
template float******** create_8d_array<float>(int len8, int len7, int len6, int len5, int len4, int len3, int len2,
                                              int len1);

template int* create_1d_array<int>(int len1);
template int** create_2d_array<int>(int len2, int len1);
template int*** create_3d_array<int>(int len3, int len2, int len1);
template int**** create_4d_array<int>(int len4, int len3, int len2, int len1);
template int***** create_5d_array<int>(int len5, int len4, int len3, int len2, int len1);
template int****** create_6d_array<int>(int len6, int len5, int len4, int len3, int len2, int len1);
template int******* create_7d_array<int>(int len7, int len6, int len5, int len4, int len3, int len2, int len1);
template int******** create_8d_array<int>(int len8, int len7, int len6, int len5, int len4, int len3, int len2,
                                          int len1);

template void delete_1d_array<float>(float *a1);
template void delete_2d_array<float>(float **a2, int len2);
template void delete_3d_array<float>(float ***a3, int len3, int len2);
template void delete_4d_array<float>(float ****a4, int len4, int len3, int len2);
template void delete_5d_array<float>(float *****a5, int len5, int len4, int len3, int len2);
template void delete_6d_array<float>(float ******a6, int len6, int len5, int len4, int len3, int len2);
template void delete_7d_array<float>(float *******a7, int len7, int len6, int len5, int len4, int len3, int len2);
template void delete_8d_array<float>(float ********a8, int len8, int len7, int len6, int len5, int len4, int len3,
                                     int len2);

template void delete_1d_array<int>(int *a1);
template void delete_2d_array<int>(int **a2, int len2);
template void delete_3d_array<int>(int ***a3, int len3, int len2);
template void delete_4d_array<int>(int ****a4, int len4, int len3, int len2);
template void delete_5d_array<int>(int *****a5, int len5, int len4, int len3, int len2);
template void delete_6d_array<int>(int ******a6, int len6, int len5, int len4, int len3, int len2);
template void delete_7d_array<int>(int *******a7, int len7, int len6, int len5, int len4, int len3, int len2);
template void delete_8d_array<int>(int ********a8, int len8, int len7, int len6, int len5, int len4, int len3,
                                   int len2);

}