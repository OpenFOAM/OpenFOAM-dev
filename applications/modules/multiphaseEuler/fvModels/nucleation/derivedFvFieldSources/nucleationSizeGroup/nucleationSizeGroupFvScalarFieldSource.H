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
    Foam::nucleationSizeGroupFvScalarFieldSource

Description
    This source condition creates a nucleation source term in the size group
    equations of a population balance. The fvSource it relates to should be a
    nucleation model.

Usage
    Example of the source condition specification:
    \verbatim
    <sourceName>
    {
        type            nucleationSizeGroup;
        libs            ("libmultiphaseEulerFvModels.so");
    }
    \endverbatim


SourceFiles
    nucleationSizeGroupFvScalarFieldSource.C

\*---------------------------------------------------------------------------*/

#ifndef nucleationSizeGroupFvScalarFieldSource_H
#define nucleationSizeGroupFvScalarFieldSource_H

#include "fvFieldSources.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class nucleationSizeGroupFvScalarFieldSource Declaration
\*---------------------------------------------------------------------------*/

class nucleationSizeGroupFvScalarFieldSource
:
    public fvScalarFieldSource
{
public:

    //- Runtime type information
    TypeName("nucleationSizeGroup");


    // Constructors

        //- Construct from internal field and dictionary
        nucleationSizeGroupFvScalarFieldSource
        (
            const DimensionedField<scalar, volMesh>&,
            const dictionary& dict
        );

        //- Copy constructor setting internal field reference
        nucleationSizeGroupFvScalarFieldSource
        (
            const nucleationSizeGroupFvScalarFieldSource&,
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
                new nucleationSizeGroupFvScalarFieldSource
                (
                    *this,
                    iF
                )
            );
        }


    //- Destructor
    virtual ~nucleationSizeGroupFvScalarFieldSource();


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
