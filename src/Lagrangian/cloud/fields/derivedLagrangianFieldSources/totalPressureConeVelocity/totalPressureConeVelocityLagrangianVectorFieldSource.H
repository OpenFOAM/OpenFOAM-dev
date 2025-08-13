/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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
    Foam::totalPressureConeVelocityLagrangianVectorFieldSource

Description
    This source condition provides a conical velocity profile, characterised by
    inner and outer cone angles, with a magnitude calculated from the
    difference between a specified upstream pressure and the local interpolated
    pressure.

Usage
    \table
        Property   | Description                          | Required? | Default
        direction  | The velocity direction               | yes       |
        thetaInner | The inner cone angle                 | yes       |
        thetaOuter | The outer cone angle                 | yes       |
        p0         | The upstream total pressure          | yes       |
        pc         | The name of the carrier pressure field | no      | p
    \endtable

    Example specification:
    \verbatim
    <LagrangianModelName>
    {
        type            totalPressureConeVelocity;
        direction       (1 0 0);
        thetaInner      15 [deg];
        thetaOuter      22.5 [deg];
        p0              1 [bar];
    }
    \endverbatim

See also
    Foam::Function1s

SourceFiles
    totalPressureConeVelocityLagrangianVectorFieldSource.C

\*---------------------------------------------------------------------------*/

#ifndef totalPressureConeVelocityLagrangianVectorFieldSource_H
#define totalPressureConeVelocityLagrangianVectorFieldSource_H

#include "LagrangianFieldSources.H"
#include "cloudLagrangianFieldSource.H"
#include "Function1LagrangianFieldSource.H"
#include "coneDirectionLagrangianVectorFieldSource.H"
#include "totalPressureVelocityMagnitudeLagrangianScalarFieldSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
      Class totalPressureConeVelocityLagrangianVectorFieldSource Declaration
\*---------------------------------------------------------------------------*/

class totalPressureConeVelocityLagrangianVectorFieldSource
:
    public LagrangianVectorFieldSource,
    private cloudLagrangianFieldSource,
    private Function1LagrangianFieldSource,
    private totalPressureVelocityMagnitudeLagrangianScalarFieldSource,
    private coneDirectionLagrangianVectorFieldSource
{
    // Private Data

        //- Direction
        const autoPtr<Function1<vector>> direction_;


    // Private Member Functions

        //- Use the base class' value method
        using Function1LagrangianFieldSource::value;


public:

    //- Runtime type information
    TypeName("totalPressureConeVelocity");


    // Constructors

        //- Construct from internal field and dictionary
        totalPressureConeVelocityLagrangianVectorFieldSource
        (
            const regIOobject&,
            const dictionary& dict
        );

        //- Copy constructor setting the internal field reference
        totalPressureConeVelocityLagrangianVectorFieldSource
        (
            const totalPressureConeVelocityLagrangianVectorFieldSource&,
            const regIOobject&
        );

        //- Construct and return a clone setting the internal field reference
        virtual autoPtr<LagrangianVectorFieldSource> clone
        (
            const regIOobject& iIo
        ) const
        {
            return autoPtr<LagrangianVectorFieldSource>
            (
                new totalPressureConeVelocityLagrangianVectorFieldSource
                (
                    *this,
                    iIo
                )
            );
        }


    //- Destructor
    virtual ~totalPressureConeVelocityLagrangianVectorFieldSource();


    // Member Functions

        //- Return the value for an instantaneous injection
        virtual tmp<LagrangianSubVectorField> value
        (
            const LagrangianInjection&,
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
