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
    Foam::Gather

Description
    Gather data from all processors onto all processors.

SourceFiles
    Gather.C

\*---------------------------------------------------------------------------*/

#ifndef Gather_H
#define Gather_H

#include "List.H"
#include "labelList.H"
#include "GatherBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class Gather Declaration
\*---------------------------------------------------------------------------*/

template<class T0>
class Gather
:
    public GatherBase,
    public List<T0>
{
    // Private data

        //- Number of processors (1 for sequential)
        label nProcs_;

        //- Storage of type 0
        // List<T0> data0_;

public:

    // Constructors

        //- Construct from containers across processors
        Gather(const T0&, const bool redistribute=true);


    // Member Functions

//        const List<T0>& data0() const
//        {
//            return data0_;
//        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "Gather.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
