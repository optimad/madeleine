/*----------------------------------------------------------------------------*\
 *
 *  G.L.O.R.I.A.
 *
 *  Optimad Engineering S.r.l. ("COMPANY") CONFIDENTIAL
 *  Copyright (c) 2014-2019 Optimad Engineering S.r.l., All Rights Reserved.
 *
 *  --------------------------------------------------------------------------
 *
 *  NOTICE:  All information contained herein is, and remains the property
 *  of COMPANY. The intellectual and technical concepts contained herein are
 *  proprietary to COMPANY and may be covered by Italian and Foreign Patents,
 *  patents in process, and are protected by trade secret or copyright law.
 *  Dissemination of this information or reproduction of this material is
 *  strictly forbidden unless prior written permission is obtained from
 *  COMPANY. Access to the source code contained herein is hereby forbidden
 *  to anyone except current COMPANY employees, managers or contractors who
 *  have executed Confidentiality and Non-disclosure agreements explicitly
 *  covering such access.
 *
 *  The copyright notice above does not evidence any actual or intended
 *  publication or disclosure of this source code, which includes information
 *  that is confidential and/or proprietary, and is a trade secret, of
 *  COMPANY. ANY REPRODUCTION, MODIFICATION, DISTRIBUTION, PUBLIC PERFORMANCE,
 *  OR PUBLIC DISPLAY OF OR THROUGH USE  OF THIS  SOURCE CODE  WITHOUT THE
 *  EXPRESS WRITTEN CONSENT OF COMPANY IS STRICTLY PROHIBITED, AND IN
 *  VIOLATION OF APPLICABLE LAWS AND INTERNATIONAL TREATIES. THE RECEIPT OR
 *  POSSESSION OF THIS SOURCE CODE AND/OR RELATED INFORMATION DOES NOT CONVEY
 *  OR IMPLY ANY RIGHTS TO REPRODUCE, DISCLOSE OR DISTRIBUTE ITS CONTENTS, OR
 *  TO MANUFACTURE, USE, OR SELL ANYTHING THAT IT  MAY DESCRIBE, IN WHOLE OR
 *  IN PART.
 *
\*----------------------------------------------------------------------------*/

#if ENABLE_MPI==1

#ifndef __MADELEINE_COMMUNICATIONS_TPP__
#define __MADELEINE_COMMUNICATIONS_TPP__

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
