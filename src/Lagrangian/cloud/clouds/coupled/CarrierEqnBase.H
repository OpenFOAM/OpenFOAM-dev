/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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
    Foam::CarrierEqnBase

Description
    Non-templated base class for carrier equations

SourceFiles
    CarrierEqnBaseI.H

\*---------------------------------------------------------------------------*/

#ifndef CarrierEqnBase_H
#define CarrierEqnBase_H

#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class CarrierEqnBase Declaration
\*---------------------------------------------------------------------------*/

class CarrierEqnBase
{
protected:

    // Protected Data

        //- Symbolic name of the equation or equation term
        word name_;

        //- Reference to the mesh
        const fvMesh& mesh_;


    // Protected Constructors

        //- Construct from components
        inline CarrierEqnBase
        (
            const word& name,
            const fvMesh& mesh
        );


public:

    // Member Functions

        //- Return the equation name
        inline const word& name() const;

        //- Return the mesh
        inline const fvMesh& mesh() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "CarrierEqnBaseI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
