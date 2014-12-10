/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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
    Foam::IDDESDelta

Description
    IDDESDelta used by the IDDES (improved low Re Spalart-Allmaras DES model)
    The min and max delta are calculated using the double distance of the min or
    max from the face centre to the cell centre.

SourceFiles
    IDDESDelta.C

\*---------------------------------------------------------------------------*/

#ifndef IDDESDeltaDelta_H
#define IDDESDeltaDelta_H

#include "incompressibleLESdelta.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class IDDESDelta Declaration
\*---------------------------------------------------------------------------*/

class IDDESDelta
:
    public LESdelta
{
    // Private data

        autoPtr<Foam::LESdelta> hmax_;
        scalar deltaCoeff_;
        scalar cw_;


    // Private Member Functions

        //- Disallow default bitwise copy construct and assignment
        IDDESDelta(const IDDESDelta&);
        void operator=(const IDDESDelta&);

        //- Calculate the delta values
        void calcDelta();


public:

    //- Runtime type information
    TypeName("IDDESDelta");


    // Constructors

        //- Construct from name, mesh and IOdictionary
        IDDESDelta
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary&
        );


    // Destructor

        ~IDDESDelta()
        {}


    // Member Functions

        //- Read the LESdelta dictionary
        void read(const dictionary&);

        // Correct values
        void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
