/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::Pstream

Description
    Inter-processor communications stream.

SourceFiles
    Pstream.C
    gatherScatter.C
    combineGatherScatter.C
    gatherScatterList.C
    exchange.C

\*---------------------------------------------------------------------------*/

#ifndef Pstream_H
#define Pstream_H

#include "UPstream.H"
#include "DynamicList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class Pstream Declaration
\*---------------------------------------------------------------------------*/

class Pstream
:
    public UPstream
{

protected:

    // Protected data

        //- Transfer buffer
        DynamicList<char> buf_;

public:

    // Declare name of the class and its debug switch
    ClassName("Pstream");


    // Constructors

        //- Construct given optional buffer size
        Pstream
        (
            const commsTypes commsType,
            const label bufSize = 0
        )
        :
            UPstream(commsType),
            buf_(0)
        {
            if (bufSize)
            {
                buf_.setCapacity(bufSize + 2*sizeof(scalar) + 1);
            }
        }


        // Gather and scatter

            //- Gather data. Apply bop to combine Value
            //  from different processors
            template<class T, class BinaryOp>
            static void gather
            (
                const List<commsStruct>& comms,
                T& Value,
                const BinaryOp& bop,
                const int tag,
                const label comm
            );

            //- Like above but switches between linear/tree communication
            template<class T, class BinaryOp>
            static void gather
            (
                T& Value,
                const BinaryOp& bop,
                const int tag = Pstream::msgType(),
                const label comm = Pstream::worldComm
            );

            //- Scatter data. Distribute without modification. Reverse of gather
            template<class T>
            static void scatter
            (
                const List<commsStruct>& comms,
                T& Value,
                const int tag,
                const label comm
            );

            //- Like above but switches between linear/tree communication
            template<class T>
            static void scatter
            (
                T& Value,
                const int tag = Pstream::msgType(),
                const label comm = Pstream::worldComm
            );

        // Combine variants. Inplace combine values from processors.
        // (Uses construct from Istream instead of <<)

            template<class T, class CombineOp>
            static void combineGather
            (
                const List<commsStruct>& comms,
                T& Value,
                const CombineOp& cop,
                const int tag,
                const label comm
            );

            //- Like above but switches between linear/tree communication
            template<class T, class CombineOp>
            static void combineGather
            (
                T& Value,
                const CombineOp& cop,
                const int tag = Pstream::msgType(),
                const label comm = Pstream::worldComm
            );

            //- Scatter data. Reverse of combineGather
            template<class T>
            static void combineScatter
            (
                const List<commsStruct>& comms,
                T& Value,
                const int tag,
                const label comm
            );

            //- Like above but switches between linear/tree communication
            template<class T>
            static void combineScatter
            (
                T& Value,
                const int tag = Pstream::msgType(),
                const label comm = Pstream::worldComm
            );

        // Combine variants working on whole List at a time.

            template<class T, class CombineOp>
            static void listCombineGather
            (
                const List<commsStruct>& comms,
                List<T>& Value,
                const CombineOp& cop,
                const int tag,
                const label comm
            );

            //- Like above but switches between linear/tree communication
            template<class T, class CombineOp>
            static void listCombineGather
            (
                List<T>& Value,
                const CombineOp& cop,
                const int tag = Pstream::msgType(),
                const label comm = Pstream::worldComm
            );

            //- Scatter data. Reverse of combineGather
            template<class T>
            static void listCombineScatter
            (
                const List<commsStruct>& comms,
                List<T>& Value,
                const int tag,
                const label comm
            );

            //- Like above but switches between linear/tree communication
            template<class T>
            static void listCombineScatter
            (
                List<T>& Value,
                const int tag = Pstream::msgType(),
                const label comm = Pstream::worldComm
            );

        // Combine variants working on whole map at a time. Container needs to
        // have iterators and find() defined.

            template<class Container, class CombineOp>
            static void mapCombineGather
            (
                const List<commsStruct>& comms,
                Container& Values,
                const CombineOp& cop,
                const int tag,
                const label comm
            );

            //- Like above but switches between linear/tree communication
            template<class Container, class CombineOp>
            static void mapCombineGather
            (
                Container& Values,
                const CombineOp& cop,
                const int tag = Pstream::msgType(),
                const label comm = UPstream::worldComm
            );

            //- Scatter data. Reverse of combineGather
            template<class Container>
            static void mapCombineScatter
            (
                const List<commsStruct>& comms,
                Container& Values,
                const int tag,
                const label comm
            );

            //- Like above but switches between linear/tree communication
            template<class Container>
            static void mapCombineScatter
            (
                Container& Values,
                const int tag = Pstream::msgType(),
                const label comm = UPstream::worldComm
            );



        // Gather/scatter keeping the individual processor data separate.
        // Values is a List of size UPstream::nProcs() where
        // Values[UPstream::myProcNo()] is the data for the current processor.

            //- Gather data but keep individual values separate
            template<class T>
            static void gatherList
            (
                const List<commsStruct>& comms,
                List<T>& Values,
                const int tag,
                const label comm
            );

            //- Like above but switches between linear/tree communication
            template<class T>
            static void gatherList
            (
                List<T>& Values,
                const int tag = Pstream::msgType(),
                const label comm = UPstream::worldComm
            );

            //- Scatter data. Reverse of gatherList
            template<class T>
            static void scatterList
            (
                const List<commsStruct>& comms,
                List<T>& Values,
                const int tag,
                const label comm
            );

            //- Like above but switches between linear/tree communication
            template<class T>
            static void scatterList
            (
                List<T>& Values,
                const int tag = Pstream::msgType(),
                const label comm = UPstream::worldComm
            );


        // Exchange

            //- Helper: exchange contiguous data. Sends sendData, receives into
            //  recvData. If block=true will wait for all transfers to finish.
            template<class Container, class T>
            static void exchange
            (
                const UList<Container>& sendData,
                const labelUList& recvSizes,
                List<Container>& recvData,
                const int tag = UPstream::msgType(),
                const label comm = UPstream::worldComm,
                const bool block = true
            );

            //- Helper: exchange sizes of sendData. sendData is the data per
            //  processor (in the communicator). Returns sizes of sendData
            //  on the sending processor.
            template<class Container>
            static void exchangeSizes
            (
                const Container& sendData,
                labelList& sizes,
                const label comm = UPstream::worldComm
            );

            //- Exchange contiguous data. Sends sendData, receives into
            //  recvData. Determines sizes to receive.
            //  If block=true will wait for all transfers to finish.
            template<class Container, class T>
            static void exchange
            (
                const UList<Container>& sendData,
                List<Container>& recvData,
                const int tag = UPstream::msgType(),
                const label comm = UPstream::worldComm,
                const bool block = true
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "gatherScatter.C"
    #include "combineGatherScatter.C"
    #include "gatherScatterList.C"
    #include "exchange.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
