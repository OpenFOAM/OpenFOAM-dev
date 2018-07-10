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
    Foam::radiation::opaqueSolid

Description
    Radiation for solid opaque solids - does nothing to energy equation source
    terms (returns zeros) but creates absorptionEmissionModel and
    scatterModel.

SourceFiles
    opaqueSolid.C

\*---------------------------------------------------------------------------*/

#ifndef opaqueSolid_H
#define opaqueSolid_H

#include "radiationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{

/*---------------------------------------------------------------------------*\
                         Class opaqueSolid Declaration
\*---------------------------------------------------------------------------*/

class opaqueSolid
:
    public radiationModel
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        opaqueSolid(const opaqueSolid&);

        //- Disallow default bitwise assignment
        void operator=(const opaqueSolid&);


public:

    //- Runtime type information
    TypeName("opaqueSolid");


    // Constructors

        //- Construct from components
        opaqueSolid(const volScalarField& T);

        //- Construct from components
        opaqueSolid(const dictionary& dict, const volScalarField& T);


    //- Destructor
    virtual ~opaqueSolid();


    // Member functions

        // Edit

            //- Solve radiation equation(s)
            void calculate();

            //- Read radiationProperties dictionary
            bool read();

            //- Source term component (for power of T^4)
            tmp<volScalarField> Rp() const;

            //- Source term component (constant)
            tmp<volScalarField::Internal> Ru() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radiation
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
