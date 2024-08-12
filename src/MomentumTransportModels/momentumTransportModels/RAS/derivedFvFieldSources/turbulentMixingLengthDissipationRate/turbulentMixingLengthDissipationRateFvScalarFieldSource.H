/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023-2024 OpenFOAM Foundation
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
    Foam::turbulentMixingLengthDissipationRateFvScalarFieldSource

Description
    This source condition provides a turbulence dissipation, \f$\epsilon\f$
    (epsilon), based on a specified mixing length. The source values are
    calculated using:

        \f[
            \epsilon = \frac{C_{\mu}^{0.75} k^{1.5}}{L}
        \f]

    where

    \vartable
        \epsilon | epsilon values
        C_{\mu}  | Model coefficient, set to 0.09
        k        | turbulence kinetic energy
        L        | length scale
    \endvartable

    In the case of a sink, the current cell values are used instead.

Usage
    \table
        Property     | Description             | Required    | Default value
        mixingLength | Length scale [m]        | yes         |
        k            | turbulence kinetic energy field name | no | k
        Cmu          | Turbulence model coefficient | no | 0.09
    \endtable

    Example of the source condition specification:
    \verbatim
    <patchName>
    {
        type            turbulentMixingLengthDissipationRate;
        mixingLength    0.005;
    }
    \endverbatim

SourceFiles
    turbulentMixingLengthDissipationRateFvScalarFieldSource.C

\*---------------------------------------------------------------------------*/

#ifndef turbulentMixingLengthDissipationRateFvScalarFieldSource_H
#define turbulentMixingLengthDissipationRateFvScalarFieldSource_H

#include "fvFieldSources.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
   Class turbulentMixingLengthDissipationRateFvScalarFieldSource Declaration
\*---------------------------------------------------------------------------*/

class turbulentMixingLengthDissipationRateFvScalarFieldSource
:
    public fvScalarFieldSource
{
private:

    // Private Data

        //- Turbulent length scale
        scalar mixingLength_;

        //- Name of the turbulent kinetic energy field
        word kName_;

        //- Turbulence model coefficient
        scalar Cmu_;


public:

    //- Runtime type information
    TypeName("turbulentMixingLengthDissipationRate");


    // Constructors

        //- Construct from internal field and dictionary
        turbulentMixingLengthDissipationRateFvScalarFieldSource
        (
            const DimensionedField<scalar, volMesh>&,
            const dictionary& dict
        );

        //- Copy constructor setting internal field reference
        turbulentMixingLengthDissipationRateFvScalarFieldSource
        (
            const turbulentMixingLengthDissipationRateFvScalarFieldSource&,
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
                new turbulentMixingLengthDissipationRateFvScalarFieldSource
                (
                    *this,
                    iF
                )
            );
        }


    //- Destructor
    virtual ~turbulentMixingLengthDissipationRateFvScalarFieldSource();


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
