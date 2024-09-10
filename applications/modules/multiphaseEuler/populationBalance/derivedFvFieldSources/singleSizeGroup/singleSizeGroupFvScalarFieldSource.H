/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024 OpenFOAM Foundation
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
    Foam::singleSizeGroupFvScalarFieldSource

Description
    Source condition for size-group fractions that sets a value of one in a
    single group, identified by its index, and a value of zero in all other
    groups

Usage
    Example of the source condition specification:
    \verbatim
    <sourceName>
    {
        type        singleSizeGroup;
        index       0;
    }
    \endverbatim

SourceFiles
    singleSizeGroupFvScalarFieldSource.C

\*---------------------------------------------------------------------------*/

#ifndef singleSizeGroupFvScalarFieldSource_H
#define singleSizeGroupFvScalarFieldSource_H

#include "fvFieldSources.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class singleSizeGroupFvScalarFieldSource Declaration
\*---------------------------------------------------------------------------*/

class singleSizeGroupFvScalarFieldSource
:
    public fvScalarFieldSource
{
private:

    // Private Data

        //- Size group index
        const scalar index_;


public:

    //- Runtime type information
    TypeName("singleSizeGroup");


    // Constructors

        //- Construct from internal field and dictionary
        singleSizeGroupFvScalarFieldSource
        (
            const DimensionedField<scalar, volMesh>&,
            const dictionary& dict
        );

        //- Copy constructor setting internal field reference
        singleSizeGroupFvScalarFieldSource
        (
            const singleSizeGroupFvScalarFieldSource&,
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
                new singleSizeGroupFvScalarFieldSource
                (
                    *this,
                    iF
                )
            );
        }


    //- Destructor
    virtual ~singleSizeGroupFvScalarFieldSource();


    // Member Functions

        //- Return the source value
        virtual tmp<scalarField> sourceValue(const fvSource&) const;

        //- Return the internal coefficient
        virtual tmp<scalarField> internalCoeff(const fvSource&) const;

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
