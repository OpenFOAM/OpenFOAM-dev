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
    Foam::coupledFacePair

Description
    Data associated with a pair of coupled faces

SourceFiles
    coupledFacePair.C

\*---------------------------------------------------------------------------*/

#ifndef coupledFacePair_H
#define coupledFacePair_H

#include "label.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

class coupledFacePair;

Ostream& operator<<(Ostream&, const coupledFacePair&);


/*---------------------------------------------------------------------------*\
                       Class coupledFacePair Declaration
\*---------------------------------------------------------------------------*/

class coupledFacePair
{
    // Private Data

        //- STAR couple ID
        label coupleIndex_;

        //- Master cell (lower cell label)
        label masterCellIndex_;

        //- Master face (lower cell label)
        label masterFaceIndex_;

        //- Slave cell (higher cell label)
        label slaveCellIndex_;

        //- Slave face (higher cell label)
        label slaveFaceIndex_;

        //- Integral or arbitrary flag
        bool integralMatch_;


public:


    // Constructors

        //- Construct from components
        coupledFacePair
        (
            const label coupleNo,              // STAR couple ID
            const label mC, const label mF,    // master cell and face
            const label sC, const label sF,    // slave cell and face
            const label integral
        );


    // Member Functions

        //- Return STAR couple ID
        label coupleIndex() const
        {
            return coupleIndex_;
        }

        //- Return master cell
        label masterCell() const
        {
            return masterCellIndex_;
        }

        //- Return master face
        label masterFace() const
        {
            return masterFaceIndex_;
        }

        //- Return slave cell
        label slaveCell() const
        {
            return slaveCellIndex_;
        }

        //- Return slave face
        label slaveFace() const
        {
            return slaveFaceIndex_;
        }

        //- Is this an integral match?
        bool integralMatch() const
        {
            return integralMatch_;
        }

    // IOstream Operators

        friend Ostream& operator<<(Ostream&, const coupledFacePair&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
