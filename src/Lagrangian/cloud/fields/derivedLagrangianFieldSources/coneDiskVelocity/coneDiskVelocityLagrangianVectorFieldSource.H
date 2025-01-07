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
    Foam::coneDiskVelocityLagrangianVectorFieldSource

Description
    This source condition provides a conical velocity profile over a disk,
    characterised by a velocity magnitude and inner and outer cone angles. The
    velocity magnitude is a Function1 and can vary with time. This condition
    must be used in conjunction with a disk injection model.

Usage
    \table
        Property      | Description                       | Required? | Default
        Umag          | The velocity magnitude            | yes       |
        thetaInner    | The inner cone angle              | yes       |
        thetaOuter    | The outer cone angle              | yes       |
    \endtable

    Example specification:
    \verbatim
    <LagrangianModelName>
    {
        type            coneDiskVelocity;
        Umag            0.4 [m/s];
        thetaInner      5 [deg];
        thetaOuter      30 [deg];
    }
    \endverbatim

SourceFiles
    coneDiskVelocityLagrangianVectorFieldSource.C

\*---------------------------------------------------------------------------*/

#ifndef coneDiskVelocityLagrangianVectorFieldSource_H
#define coneDiskVelocityLagrangianVectorFieldSource_H

#include "LagrangianFieldSources.H"
#include "Function1LagrangianFieldSource.H"
#include "randomGenerator.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class coneDiskVelocityLagrangianVectorFieldSource Declaration
\*---------------------------------------------------------------------------*/

class coneDiskVelocityLagrangianVectorFieldSource
:
    public LagrangianVectorFieldSource,
    private Function1LagrangianFieldSource<vector>
{
private:

    // Private Data

        //- Velocity magnitude
        const autoPtr<Function1<scalar>> Umag_;

        //- Inner half-cone angle
        const autoPtr<Function1<scalar>> thetaInner_;

        //- Outer half-cone angle
        const autoPtr<Function1<scalar>> thetaOuter_;


    // Private Member Functions

        //- Use the base class' value method
        using Function1LagrangianFieldSource<vector>::value;


public:

    //- Runtime type information
    TypeName("coneDiskVelocity");


    // Constructors

        //- Construct from internal field and dictionary
        coneDiskVelocityLagrangianVectorFieldSource
        (
            const regIOobject&,
            const dictionary& dict
        );

        //- Copy constructor setting the internal field reference
        coneDiskVelocityLagrangianVectorFieldSource
        (
            const coneDiskVelocityLagrangianVectorFieldSource&,
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
                new coneDiskVelocityLagrangianVectorFieldSource(*this, iIo)
            );
        }


    //- Destructor
    virtual ~coneDiskVelocityLagrangianVectorFieldSource();


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
