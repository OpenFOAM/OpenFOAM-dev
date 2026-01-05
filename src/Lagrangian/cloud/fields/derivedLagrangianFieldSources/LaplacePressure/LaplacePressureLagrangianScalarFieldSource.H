/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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
    Foam::LaplacePressureLagrangianScalarFieldSource

Description
    This source condition sets the values of the pressure to that of the
    carrier plus an amount caused by the effect of the surface tension on a
    spherical droplet or bubble.

    Note that for now, the surface tension is an input parameter. This is
    likely to change as more surface-tension-dependent models are introduced
    and it becomes necessary to centralise the specification of the surface
    tension modelling. For now, this serves as a proof of concept that the
    pressure of the particles can differ from that of the carrier.

Usage
    \table
        Property | Description                     | Required? | Default
        sigma    | The surface tension coefficient | yes       |
    \endtable

    Example specification:
    \verbatim
    <LagrangianModelName>
    {
        type            LaplacePressure;
        sigma           0.07 [N/m];
    }
    \endverbatim

SourceFiles
    LaplacePressureLagrangianScalarFieldSource.C

\*---------------------------------------------------------------------------*/

#ifndef LaplacePressureLagrangianScalarFieldSource_H
#define LaplacePressureLagrangianScalarFieldSource_H

#include "LagrangianFieldSources.H"
#include "cloudLagrangianFieldSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
       Class LaplacePressureLagrangianScalarFieldSource Declaration
\*---------------------------------------------------------------------------*/

class LaplacePressureLagrangianScalarFieldSource
:
    public LagrangianScalarFieldSource,
    private cloudLagrangianFieldSource
{
private:

    // Private Data

        //- The surface tension coefficient
        const dimensionedScalar sigma_;


public:

    //- Runtime type information
    TypeName("LaplacePressure");


    // Constructors

        //- Construct from internal field and dictionary
        LaplacePressureLagrangianScalarFieldSource
        (
            const regIOobject&,
            const dictionary& dict
        );

        //- Copy constructor setting the internal field reference
        LaplacePressureLagrangianScalarFieldSource
        (
            const LaplacePressureLagrangianScalarFieldSource&,
            const regIOobject&
        );

        //- Construct and return a clone setting the internal field reference
        virtual autoPtr<LagrangianScalarFieldSource> clone
        (
            const regIOobject& iIo
        ) const
        {
            return autoPtr<LagrangianScalarFieldSource>
            (
                new LaplacePressureLagrangianScalarFieldSource(*this, iIo)
            );
        }


    //- Destructor
    virtual ~LaplacePressureLagrangianScalarFieldSource();


    // Member Functions

        //- Return the value for an instantaneous injection
        virtual tmp<LagrangianSubScalarField> value
        (
            const LagrangianInjection&,
            const LagrangianSubMesh&
        ) const;

        //- Return the value without a model
        virtual tmp<LagrangianSubScalarField> value
        (
            const LagrangianSubMesh&
        ) const;

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
