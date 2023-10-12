/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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
    Foam::energyFvScalarFieldSource

Description
    This source condition provides a value for the energy derived from the
    corresponding condition for the temperature, and the conditions for other
    relevant thermodynamic variables. This is constructed automatically by the
    thermodynamic model. The user does not specify it.

SourceFiles
    energyFvScalarFieldSource.C

\*---------------------------------------------------------------------------*/

#ifndef energyFvScalarFieldSource_H
#define energyFvScalarFieldSource_H

#include "fvFieldSources.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class energyFvScalarFieldSource Declaration
\*---------------------------------------------------------------------------*/

class energyFvScalarFieldSource
:
    public fvScalarFieldSource
{
public:

    //- Runtime type information
    TypeName("energy");


    // Constructors

        //- Construct from internal field
        energyFvScalarFieldSource(const DimensionedField<scalar, volMesh>&);

        //- Construct from internal field and dictionary
        energyFvScalarFieldSource
        (
            const DimensionedField<scalar, volMesh>&,
            const dictionary& dict
        );

        //- Copy constructor setting internal field reference
        energyFvScalarFieldSource
        (
            const energyFvScalarFieldSource&,
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
                new energyFvScalarFieldSource(*this, iF)
            );
        }


    //- Destructor
    virtual ~energyFvScalarFieldSource();


    // Member Functions

        //- Return the source value
        virtual tmp<scalarField> sourceValue(const fvSource&) const;

        //- Return the internal coefficient
        virtual tmp<scalarField> internalCoeff(const fvSource&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
