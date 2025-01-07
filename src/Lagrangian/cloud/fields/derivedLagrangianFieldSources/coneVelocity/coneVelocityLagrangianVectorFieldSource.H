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
    Foam::coneVelocityLagrangianVectorFieldSource

Description
    This source condition provides a conical velocity profile, characterised by
    a mean velocity and inner and outer cone angles. The mean velocity is a
    Function1 and can vary with time. This is likely to be used in conjunction
    with a point injection model.

Usage
    \table
        Property      | Description                       | Required? | Default
        Umean         | The velocity along the centreline | yes       |
        thetaInner    | The inner cone angle              | yes       |
        thetaOuter    | The outer cone angle              | yes       |
    \endtable

    Example specification:
    \verbatim
    <LagrangianModelName>
    {
        type            coneVelocity;
        Umean
        {
            type            table;
            values
            (
                (0 (0.1 0 0))
                (0.5 (0.5 0 0))
            );
        }
        thetaInner      15 [deg];
        thetaOuter      22.5 [deg];
    }
    \endverbatim

See also
    Foam::Function1s

SourceFiles
    coneVelocityLagrangianVectorFieldSource.C

\*---------------------------------------------------------------------------*/

#ifndef coneVelocityLagrangianVectorFieldSource_H
#define coneVelocityLagrangianVectorFieldSource_H

#include "LagrangianFieldSources.H"
#include "Function1LagrangianFieldSource.H"
#include "restartableRandomGenerator.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class coneVelocityLagrangianVectorFieldSource Declaration
\*---------------------------------------------------------------------------*/

class coneVelocityLagrangianVectorFieldSource
:
    public LagrangianVectorFieldSource,
    private Function1LagrangianFieldSource<vector>
{
private:

    // Private Data

        //- Mean velocity
        const autoPtr<Function1<vector>> Umean_;

        //- Inner half-cone angle
        const autoPtr<Function1<scalar>> thetaInner_;

        //- Outer half-cone angle
        const autoPtr<Function1<scalar>> thetaOuter_;

        //- A random number generator
        mutable restartableRandomGenerator rndGen_;

        //- The time index
        mutable label timeIndex_;


    // Private Member Functions

        //- Use the base class' value method
        using Function1LagrangianFieldSource<vector>::value;


public:

    //- Runtime type information
    TypeName("coneVelocity");


    // Constructors

        //- Construct from internal field and dictionary
        coneVelocityLagrangianVectorFieldSource
        (
            const regIOobject&,
            const dictionary& dict
        );

        //- Copy constructor setting the internal field reference
        coneVelocityLagrangianVectorFieldSource
        (
            const coneVelocityLagrangianVectorFieldSource&,
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
                new coneVelocityLagrangianVectorFieldSource(*this, iIo)
            );
        }


    //- Destructor
    virtual ~coneVelocityLagrangianVectorFieldSource();


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
