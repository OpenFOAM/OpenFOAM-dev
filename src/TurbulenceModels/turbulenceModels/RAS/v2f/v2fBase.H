/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2018 OpenFOAM Foundation
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
    Foam::RASModels::v2fBase

Description
    Abstract base-class for v2-f models to provide BCs access to the v2 and f
    fields.

See also
    Foam::RASModels::v2f

SourceFiles
    v2fBase.C

\*---------------------------------------------------------------------------*/

#ifndef v2fBase_H
#define v2fBase_H

#include "RASModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                            Class v2fBase Declaration
\*---------------------------------------------------------------------------*/

class v2fBase
{
public:

    //- Runtime type information
    TypeName("v2fBase");


    // Constructors

        v2fBase()
        {}


    //- Destructor
    virtual ~v2fBase()
    {}


    // Member Functions

        //- Return turbulence stress normal to streamlines
        virtual tmp<volScalarField> v2() const = 0;

        //- Return the damping function
        virtual tmp<volScalarField> f() const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
