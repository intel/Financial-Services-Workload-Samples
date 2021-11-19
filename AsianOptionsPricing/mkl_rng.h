#include <mkl.h>
#include <cassert>
#include <memory>

template <typename T = double>
class mkl_rng
{
public:
    mkl_rng(const int brng, 
        const unsigned int seed, 
        const int _distribution, 
        const T _mean,
        const T _sigma,
        const size_t _buffer_size,
        const size_t _aligment
    ) :
        distribution(_distribution), mean(_mean), sigma(_sigma), buffer_size(_buffer_size), aligment(_aligment), state_ptr(NULL), data(NULL)
    {  
        if (vslNewStream(&state_ptr, brng, seed) == VSL_STATUS_OK)
        {
#if 0
            auto raw_buffer = buffer_size + aligment;
            raw_data = malloc(raw_buffer);
            auto raw_data_tmp = raw_data;
            data = (T*)std::align(aligment, buffer_size, raw_data, raw_buffer);
            raw_data = raw_data_tmp;
#endif
            data = (T*)mkl_malloc(buffer_size, aligment);
        }
            
    }

    mkl_rng(const mkl_rng& _rng) : 
        state_ptr(_rng.state_ptr), distribution(_rng.distribution), mean(_rng.mean), sigma(_rng.sigma), buffer_size(_rng.buffer_size), aligment(_rng.aligment)
    {
        data = (T*) mkl_malloc(buffer_size, aligment);
    }
    
    T* get_gaussian(int);

    ~mkl_rng()
    {
        mkl_free(data);
    }

private:
    mkl_rng() {};
    VSLStreamStatePtr state_ptr;
    int distribution;
    const size_t buffer_size;
    const size_t aligment;
    T mean;
    T sigma;
    T* data;
};

inline double* mkl_rng<double>::get_gaussian(int rngNums)
{
    assert(buffer_size >= sizeof(double)* rngNums);
    vdRngGaussian(distribution, state_ptr, rngNums, data, mean, sigma);
    return data;
}

inline float* mkl_rng<float>::get_gaussian(int rngNums)
{
    assert(buffer_size >= sizeof(float) * rngNums);
    vsRngGaussian(distribution, state_ptr, rngNums, data, mean, sigma);
    return data;
}
