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
    Foam::singleGroupFractionFvScalarFieldSource

Description
    Source condition for group fractions of a population balance model that
    sets a value of one in a single group, identified by its index, and a value
    of zero in all other groups

Usage
    Example of the source condition specification:
    \verbatim
    <sourceName>
    {
        type        singleGroupFraction;
        index       0;
    }
    \endverbatim

SourceFiles
    singleGroupFractionFvScalarFieldSource.C

\*---------------------------------------------------------------------------*/

#ifndef singleGroupFractionFvScalarFieldSource_H
#define singleGroupFractionFvScalarFieldSource_H

#include "fvFieldSources.H"
#include "groupPropertyFvScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class singleGroupFractionFvScalarFieldSource Declaration
\*---------------------------------------------------------------------------*/

class singleGroupFractionFvScalarFieldSource
:
    public fvScalarFieldSource,
    public groupPropertyFvScalarField
{
private:

    // Private Data

        //- Group index
        const scalar index_;


    // Private Member Functions

        //- Return the allocation coefficient value
        scalar eta() const;


public:

    //- Runtime type information
    TypeName("singleGroupFraction");


    // Constructors

        //- Construct from internal field and dictionary
        singleGroupFractionFvScalarFieldSource
        (
            const DimensionedField<scalar, volMesh>&,
            const dictionary& dict
        );

        //- Copy constructor setting internal field reference
        singleGroupFractionFvScalarFieldSource
        (
            const singleGroupFractionFvScalarFieldSource&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual autoPtr<fvScalarFieldSource> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return autoPtr<fvScalarFieldSource>
            (
                new singleGroupFractionFvScalarFieldSource
                (
                    *this,
                    iF
                )
            );
        }


    //- Destructor
    virtual ~singleGroupFractionFvScalarFieldSource();


    // Member Functions

        //- Return the source value
        virtual tmp<DimensionedField<scalar, volMesh>> sourceValue
        (
            const fvSource& model,
            const DimensionedField<scalar, volMesh>& source
        ) const;

        //- Return the source value
        virtual tmp<scalarField> sourceValue
        (
            const fvSource& model,
            const scalarField& source,
            const labelUList& cells
        ) const;

        //- Return the internal coefficient
        virtual tmp<DimensionedField<scalar, volMesh>> internalCoeff
        (
            const fvSource& model,
            const DimensionedField<scalar, volMesh>& source
        ) const;

        //- Return the internal coefficient
        virtual tmp<scalarField> internalCoeff
        (
            const fvSource& model,
            const scalarField& source,
            const labelUList& cells
        ) const;

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
