/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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
    Foam::polyRemovePoint

Description
    Class containing data for point removal.

\*---------------------------------------------------------------------------*/

#ifndef polyRemovePoint_H
#define polyRemovePoint_H

#include "label.H"
#include "topoAction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class polyRemovePoint Declaration
\*---------------------------------------------------------------------------*/

class polyRemovePoint
:
    public topoAction
{
    // Private Data

        //- Point ID
        label pointIndex_;

        //- Merge point ID or -1
        label mergePointIndex_;


public:

    // Static Data Members

        //- Runtime type information
        TypeName("removePoint");


    // Constructors

        //- Construct null.  Used for constructing lists
        polyRemovePoint()
        :
            pointIndex_(-1),
            mergePointIndex_(-1)
        {}

        //- Construct from components
        polyRemovePoint(const label pointID, const label mergePointID = -1)
        :
            pointIndex_(pointID),
            mergePointIndex_(mergePointID)
        {}

        //- Construct and return a clone
        virtual autoPtr<topoAction> clone() const
        {
            return autoPtr<topoAction>(new polyRemovePoint(*this));
        }


    // Default Destructor


    // Member Functions

        //- Return point ID
        label pointIndex() const
        {
            return pointIndex_;
        }

        label mergePointIndex() const
        {
            return mergePointIndex_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
