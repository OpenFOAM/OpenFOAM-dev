/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019-2020 OpenFOAM Foundation
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
    Foam::dragModels::aerosolDrag

Description
    This drag model is intended for simulations involving dispersed phase
    with a very small diameter. It combines the Stokes drag with a
    Cunningham correction factor.

SourceFiles
    aerosolDrag.C

\*---------------------------------------------------------------------------*/

#ifndef aerosolDrag_H
#define aerosolDrag_H

#include "dragModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace dragModels
{

/*---------------------------------------------------------------------------*\
                           Class aerosolDrag Declaration
\*---------------------------------------------------------------------------*/

class aerosolDrag
:
    public dragModel
{
    // Private Data

        //- Model coefficient A1
        const scalar A1_;

        //- Model coefficient A2
        const scalar A2_;

        //- Model coefficient A3
        const scalar A3_;

        //- Diameter of the continuous phase molecules, defaults to value for N2
        const dimensionedScalar dm_;


public:

    //- Runtime type information
    TypeName("aerosolDrag");


    // Constructors

        //- Construct from a dictionary and a phase pair
        aerosolDrag
        (
            const dictionary& dict,
            const phasePair& pair,
            const bool registerObject
        );


    //- Destructor
    virtual ~aerosolDrag();


    // Member Functions

        //- Drag coefficient
        virtual tmp<volScalarField> CdRe() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace dragModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
