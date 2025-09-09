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
    Foam::secondaryPropertyFvScalarFieldSource

Description
    Base class for source conditions of secondary properties of the size-groups
    of a population balance model.

SourceFiles
    secondaryPropertyFvScalarFieldSource.C

\*---------------------------------------------------------------------------*/

#ifndef secondaryPropertyFvScalarFieldSource_H
#define secondaryPropertyFvScalarFieldSource_H

#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
namespace diameterModels
{
    class sizeGroup;
}

/*---------------------------------------------------------------------------*\
          Class secondaryPropertyFvScalarFieldSource Declaration
\*---------------------------------------------------------------------------*/

class secondaryPropertyFvScalarFieldSource
{
private:

    // Private Data

        //- Reference to the internal field
        const DimensionedField<scalar, volMesh>& internalField_;

        //- The index of this group
        mutable label i_;


public:

    // Constructors

        //- Construct from internal field
        secondaryPropertyFvScalarFieldSource
        (
            const DimensionedField<scalar, volMesh>&
        );


public:

    // Member Functions

        //- Return the list of size-groups
        const UPtrList<diameterModels::sizeGroup>& fis() const;

        //- Return the index of this group
        label i() const;

        //- Return a size-group
        const diameterModels::sizeGroup& fi(const label deltai = 0) const;

        //- Return a secondary property field
        template<class ModelType>
        const volScalarField& fld(const label deltai = 0) const;

        //- Return a secondary property model
        template<class ModelType>
        const ModelType& model(const label deltai = 0) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "secondaryPropertyFvScalarFieldSourceTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
