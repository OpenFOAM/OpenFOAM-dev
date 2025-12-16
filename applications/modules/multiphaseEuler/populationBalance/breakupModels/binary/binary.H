/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2025 OpenFOAM Foundation
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
    Foam::populationBalance::breakupModels::binary

Description
    Base class for binary breakup models that provide a breakup rate between a
    size class pair directly, i.e. without explicitly stating the daughter size
    distribution function.

See also
    Foam::populationBalanceModel

SourceFiles
    binary.C

\*---------------------------------------------------------------------------*/

#ifndef binary_H
#define binary_H

#include "breakupModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalance
{
namespace breakupModels
{

/*---------------------------------------------------------------------------*\
                     Class binary Declaration
\*---------------------------------------------------------------------------*/

class binary
:
    public breakupModel
{
public:

    // Constructor

        //- Inherit constructors
        using breakupModel::breakupModel;


    //- Destructor
    virtual ~binary()
    {}


    // Member Functions

        //- Return the breakup rate between two groups
        virtual tmp<volScalarField::Internal> rate
        (
            const label i,
            const label j
        ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace breakupModels
} // End namespace populationBalance
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
