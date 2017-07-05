#ifndef MIKAN_MK_MEMORY_HPP_
#define MIKAN_MK_MEMORY_HPP_

namespace mikan {

template<typename T>
T *create_1d_array(int len1);

template<typename T>
T **create_2d_array(int len2, int len1);

template<typename T>
T ***create_3d_array(int len3, int len2, int len1);

template<typename T>
T ****create_4d_array(int len4, int len3, int len2, int len1);

template<typename T>
T *****create_5d_array(int len5, int len4, int len3, int len2, int len1);

template<typename T>
T ******create_6d_array(int len6, int len5, int len4, int len3, int len2, int len1);

template<typename T>
T *******create_7d_array(int len7, int len6, int len5, int len4, int len3, int len2, int len1);

template<typename T>
T ********create_8d_array(int len8, int len7, int len6, int len5, int len4, int len3, int len2, int len1);

template<typename T>
void delete_1d_array(T *a1);

template<typename T>
void delete_2d_array(T **a2, int len2);

template<typename T>
void delete_3d_array(T ***a3, int len3, int len2);

template<typename T>
void delete_4d_array(T ****a4, int len4, int len3, int len2);

template<typename T>
void delete_5d_array(T *****a5, int len5, int len4, int len3, int len2);

template<typename T>
void delete_6d_array(T ******a6, int len6, int len5, int len4, int len3, int len2);

template<typename T>
void delete_7d_array(T *******a7, int len7, int len6, int len5, int len4, int len3, int len2);

template<typename T>
void delete_8d_array(T ********a8, int len8, int len7, int len6, int len5, int len4, int len3, int len2);

}

#endif //MIKAN_MK_MEMORY_HPP_
