#pragma once

namespace ctrs
{
    template<typename dtype, const size_t ar_size> struct array
    {
        dtype data[ar_size];
        dtype& operator [] (size_t idx) {return data[idx];}
        dtype* begin() noexcept {return &data[0];}
        dtype* end()   noexcept {return &data[0]+ar_size;}
        constexpr size_t size(void) const noexcept {return ar_size;}
        void fill(const dtype& val)
        {
            for (size_t i = 0; i < this->size(); i++) data[i] = val;
        }
        array(const dtype& val) {fill(val);}
        array(void) {fill(dtype());}
    };
}