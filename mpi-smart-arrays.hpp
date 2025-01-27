#ifndef _MPI_SMART_ARRAYS_HPP_
#define _MPI_SMART_ARRAYS_HPP_

#include <exception>
#include <iostream>
#include <mpi.h>
#include <span>
#include <vector>

namespace smartarray {

template <typename T> class SmartArray {
  private:
    std::size_t m_global_size;
    std::size_t m_elem_size;
    MPI_Comm m_comm;
    std::size_t m_max_window_size;
    int m_comm_rank;
    int m_comm_size;
    int m_right_neighbor;
    int m_left_neighbor;
    std::size_t m_local_size;
    std::size_t m_my_local_size;
    std::vector<T> m_local_buffer;

  public:
    SmartArray(std::size_t global_size, std::size_t elem_size, MPI_Comm comm,
               std::size_t max_window_size = 10)
        : m_global_size(global_size), m_elem_size(elem_size), m_comm(comm),
          m_max_window_size(max_window_size) {
        MPI_Comm_rank(m_comm, &m_comm_rank);
        MPI_Comm_size(m_comm, &m_comm_size);
        if (static_cast<std::size_t>(m_comm_size) > global_size) {
            throw std::runtime_error("Array is smaller than number of ranks");
        }
        if (m_elem_size == 0) {
            throw std::runtime_error("elem_size can't be smaller than 1");
        }
        m_right_neighbor = (m_comm_rank + 1) % m_comm_size;
        m_left_neighbor = ((m_comm_rank - 1) + m_comm_size) % m_comm_size;
        m_local_size = (m_global_size + (m_comm_size - 1)) / m_comm_size;

        m_max_window_size = std::min(m_local_size, m_max_window_size);

        m_my_local_size = m_local_size;

        if (m_comm_rank == m_comm_size - 1) {
            m_my_local_size = m_global_size - m_local_size * (m_comm_size - 1);
        }

        m_local_buffer.resize((m_my_local_size + 2 * m_max_window_size) *
                              m_elem_size);
    }

    void print() {
        for (int i = 0; i < m_comm_size; i++) {
            if (i == m_comm_rank) {
                for (int j = 0; j < m_my_local_size; j++) {
                    auto* buffer_start = m_local_buffer.data() +
                                         (m_max_window_size + j) * m_elem_size;
                    for (int k = 0; k < m_elem_size; k++) {
                        std::cout << buffer_start[k] << " ";
                    }
                    std::cout << std::endl;
                }
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }
    }

    void sync() {
        std::size_t right_send_start_idx = m_my_local_size * m_elem_size;
        std::size_t left_recv_start_idx = 0;

        MPI_Sendrecv(m_local_buffer.data() + right_send_start_idx,
                     sizeof(T) * m_max_window_size * m_elem_size, MPI_BYTE,
                     m_right_neighbor, 0,
                     m_local_buffer.data() + left_recv_start_idx,
                     sizeof(T) * m_max_window_size * m_elem_size, MPI_BYTE,
                     m_left_neighbor, 0, m_comm, MPI_STATUS_IGNORE);

        std::size_t left_send_start_idx = m_max_window_size * m_elem_size;
        std::size_t right_recv_start_idx =
            (m_my_local_size + m_max_window_size) * m_elem_size;

        MPI_Sendrecv(m_local_buffer.data() + left_send_start_idx,
                     sizeof(T) * m_max_window_size * m_elem_size, MPI_BYTE,
                     m_left_neighbor, 0,
                     m_local_buffer.data() + right_recv_start_idx,
                     sizeof(T) * m_max_window_size * m_elem_size, MPI_BYTE,
                     m_right_neighbor, 0, m_comm, MPI_STATUS_IGNORE);
    }

    template <typename Transform> void transform(Transform f) {
#pragma omp parallel for
        for (std::size_t i = 0; i < m_my_local_size; i++) {
            f(m_local_buffer.data() + (i + m_max_window_size) * m_elem_size,
              m_comm_rank * m_local_size + i, m_global_size, m_elem_size);
        }
    }

    template <typename Transform>
    void transform_window(Transform f, std::size_t window_size) {
        if (!(window_size <= m_max_window_size)) {
            throw std::runtime_error(
                "window_size=" + std::to_string(window_size) +
                " is > max_window_size=" + std::to_string(m_max_window_size));
        }
        sync();
#pragma omp parallel for
        for (std::size_t i = 0; i < m_my_local_size; i++) {
            const auto* window_start =
                m_local_buffer.data() + (m_max_window_size + i) * m_elem_size;
            f(m_local_buffer.data() + (i + m_max_window_size) * m_elem_size,
              window_start, window_size, m_comm_rank * m_local_size + i,
              m_global_size, m_elem_size);
        }
    }
};

} // namespace smartarray

#endif // _MPI_SMART_ARRAYS_HPP_