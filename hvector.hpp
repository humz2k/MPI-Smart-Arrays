#include <stdio.h>
#include <stdlib.h>

template<class T>
struct hvec3{
    T x;
    T y;
    T z;

    inline hvec3<T> operator+(const hvec3<T> other){
        hvec3<T> out;
        out.x = x + other.x;
        out.y = y + other.y;
        out.z = z + other.z;
        return out;
    }

    inline hvec3<T> operator+(const T other){
        hvec3<T> out;
        out.x = x + other;
        out.y = y + other;
        out.z = z + other;
        return out;
    }

    inline hvec3<T> operator-(const hvec3<T> other){
        hvec3<T> out;
        out.x = x - other.x;
        out.y = y - other.y;
        out.z = z - other.z;
        return out;
    }

    inline hvec3<T> operator-(const T other){
        hvec3<T> out;
        out.x = x - other;
        out.y = y - other;
        out.z = z - other;
        return out;
    }

    inline hvec3<T> operator*(const hvec3<T> other){
        hvec3<T> out;
        out.x = x * other.x;
        out.y = y * other.y;
        out.z = z * other.z;
        return out;
    }

    inline hvec3<T> operator*(const T other){
        hvec3<T> out;
        out.x = x * other;
        out.y = y * other;
        out.z = z * other;
        return out;
    }

    inline hvec3<T> operator/(const hvec3<T> other){
        hvec3<T> out;
        out.x = x / other.x;
        out.y = y / other.y;
        out.z = z / other.z;
        return out;
    }

    inline hvec3<T> operator/(const T other){
        hvec3<T> out;
        out.x = x / other;
        out.y = y / other;
        out.z = z / other;
        return out;
    }

    T& operator[](int idx){
        if (idx == 0){
            return x;
        }
        if (idx == 1){
            return y;
        }
        return z;
    }
};

template<class T>
inline hvec3<T> make_hvec3(T x, T y, T z){
    hvec3<T> out;
    out.x = x;
    out.y = y;
    out.z = z;
    return out;
}

