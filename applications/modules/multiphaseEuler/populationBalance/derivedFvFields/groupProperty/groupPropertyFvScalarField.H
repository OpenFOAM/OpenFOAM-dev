/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024-2025 OpenFOAM Foundation
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
    Foam::groupPropertyFvScalarField

Description
    Base class for source conditions of properties of the groups in a
    population balance model.

SourceFiles
    groupPropertyFvScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef groupPropertyFvScalarField_H
#define groupPropertyFvScalarField_H

#include "DimensionedField.H"
#include "volMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class populationBalanceModel;

/*---------------------------------------------------------------------------*\
                Class groupPropertyFvScalarField Declaration
\*---------------------------------------------------------------------------*/

class groupPropertyFvScalarField
{
private:

    // Private Data

        //- Reference to the internal field
        const DimensionedField<scalar, volMesh>& internalField_;

        //- Pointer to the population balance
        mutable const populationBalanceModel* popBalPtr_;

        //- The index of this group
        mutable label i_;


public:

    // Constructors

        //- Construct from internal field
        groupPropertyFvScalarField
        (
            const DimensionedField<scalar, volMesh>&
        );


public:

    // Member Functions

        //- Return the population balance
        const populationBalanceModel& popBal() const;

        //- Return the index of this group
        label i() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
