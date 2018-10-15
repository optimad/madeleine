#if ENABLE_MPI==1

#ifndef __GLORIA_COMMUNICATIONS_TPP__
#define __GLORIA_COMMUNICATIONS_TPP__

#include <bitpit_IO.hpp>

/*!
    \class ListBufferStreamer

    \brief The ListBufferStreamer class allows to stream list data from / to
    the buffer of a ListCommunicator.
*/

/*!
    Creates a new streamer

    \param container is the container that holds the data that will be
    exchanged
*/
template<typename container_t, typename value_t>
ListBufferStreamer<container_t, value_t>::ListBufferStreamer(container_t *container)
    : ExchangeBufferStreamer(sizeof(value_t)),
      m_container(container)
{
}

/*!
    Creates a new streamer

    \param container is the container that holds the data that will be
    exchanged
    \param itemSize is the size, expressed in bytes, of the single item that
    will be exchanged
*/
template<typename container_t, typename value_t>
ListBufferStreamer<container_t, value_t>::ListBufferStreamer(container_t *container, const size_t &itemSize)
    : ExchangeBufferStreamer(itemSize),
      m_container(container)
{
}

/*!
    Read the dataset from the buffer.

    \param rank is the rank of the processor who sent the data
    \param buffer is the buffer where the data will be read from
    \param list is the list of ids that will be read
*/
template<typename container_t, typename value_t>
void ListBufferStreamer<container_t, value_t>::read(int const &rank, bitpit::RecvBuffer &buffer,
                                      const std::vector<long> &list)
{
    BITPIT_UNUSED(rank);

    for (const long k : list) {
        buffer >> (*m_container)[k];
    }
}

/*!
    Write the dataset to the buffer.

    \param rank is the rank of the processor who will receive the data
    \param buffer is the buffer where the data will be written to
    \param list is the list of ids that will be written
*/
template<typename container_t, typename value_t>
void ListBufferStreamer<container_t, value_t>::write(const int &rank, bitpit::SendBuffer &buffer,
                                                    const std::vector<long> &list)
{
    BITPIT_UNUSED(rank);

    for (const long k : list) {
        buffer << (*m_container)[k];
    }
}

/*!
    Gets a reference to the container that holds the data that will be
    streamed.

    \result A reference to the container that holds the data that will be
    streamed.
*/
template<typename container_t, typename value_t>
container_t & ListBufferStreamer<container_t, value_t>::getContainer()
{
    return *m_container;
}

#endif

#endif
