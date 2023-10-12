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
    Foam::turbulentIntensityKineticEnergyFvScalarFieldSource

Description
    This source condition provides a turbulent kinetic energy based on
    user-supplied turbulence intensity, defined as a fraction of the
    mean velocity:

        \f[
            k = 1.5 (I |U|)^2
        \f]

    where

    \vartable
        k       | kinetic energy values
        I       | turbulence intensity
        U       | velocity field
    \endvartable

    In the case of a sink, the current cell values are used instead.

Usage
    \table
        Property     | Description             | Required    | Default value
        intensity    | fraction of mean field [0-1] | yes    |
        U            | velocity field name     | no          | U
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type        turbulentIntensityKineticEnergyInlet;
        intensity   0.05;           // 5% turbulence
    }
    \endverbatim

SourceFiles
    turbulentIntensityKineticEnergyFvScalarFieldSource.C

\*---------------------------------------------------------------------------*/

#ifndef turbulentIntensityKineticEnergyFvScalarFieldSource_H
#define turbulentIntensityKineticEnergyFvScalarFieldSource_H

#include "fvFieldSources.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
     Class turbulentIntensityKineticEnergyFvScalarFieldSource Declaration
\*---------------------------------------------------------------------------*/

class turbulentIntensityKineticEnergyFvScalarFieldSource
:
    public fvScalarFieldSource
{
private:

    // Private Data

        //- Turbulent intensity as fraction of mean velocity
        scalar intensity_;

        //- Name of the velocity field
        word UName_;


public:

    //- Runtime type information
    TypeName("turbulentIntensityKineticEnergy");


    // Constructors

        //- Construct from internal field and dictionary
        turbulentIntensityKineticEnergyFvScalarFieldSource
        (
            const DimensionedField<scalar, volMesh>&,
            const dictionary& dict
        );

        //- Copy constructor setting internal field reference
        turbulentIntensityKineticEnergyFvScalarFieldSource
        (
            const turbulentIntensityKineticEnergyFvScalarFieldSource&,
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
                new turbulentIntensityKineticEnergyFvScalarFieldSource
                (
                    *this,
                    iF
                )
            );
        }


    //- Destructor
    virtual ~turbulentIntensityKineticEnergyFvScalarFieldSource();


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
