#ifndef _MPI_SMART_ARRAYS_HPP_
#define _MPI_SMART_ARRAYS_HPP_

#include <vector>
#include <mpi.h>
#include <iostream>

namespace smartarray {

template<typename T>
class SmartArray {
private:
    std::size_t m_global_size;
    MPI_Comm m_comm;
    int m_comm_rank;
    int m_comm_size;
    std::size_t m_local_size;
    std::vector<T> m_local_buffer;
public:
    SmartArray(std::size_t global_size, MPI_Comm comm) : m_global_size(global_size), m_comm(comm){
        MPI_Comm_rank(m_comm, &m_comm_rank);
        MPI_Comm_size(m_comm, &m_comm_size);
        m_local_size = (m_global_size + (m_comm_size - 1)) / m_comm_size;
        if (!m_comm_rank){
            std::cout << "local_size = " << m_local_size << std::endl;
        }
    }
};

} // namespace smartarray

#endif // _MPI_SMART_ARRAYS_HPP_