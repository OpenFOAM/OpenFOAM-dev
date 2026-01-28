/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025-2026 OpenFOAM Foundation
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
    Foam::fanVelocityLagrangianVectorFieldSource

Description
    This source condition provides a fan velocity profile, characterised by a
    centreline velocity, a normal, and inner and outer fan angles. The
    centreline velocity and normal are Function1s and can vary with time. This
    is likely to be used in conjunction with a point injection model.

Usage
    \table
        Property      | Description                       | Required? | Default
        Ucentre       | The velocity along the centreline | yes       |
        normal        | Normal to the plane of the fan    | yes       |
        thetaInner    | The inner fan angle               | yes       |
        thetaOuter    | The outer fan angle               | yes       |
    \endtable

    Example specification:
    \verbatim
    <LagrangianModelName>
    {
        type            fanVelocity;
        Ucentre
        {
            type            table;
            values
            (
                (0 (0.1 0 0))
                (0.5 (0.5 0 0))
            );
        }
        normal          (0 0 1);
        thetaInner      15 [deg];
        thetaOuter      22.5 [deg];
    }
    \endverbatim

See also
    Foam::Function1s

SourceFiles
    fanVelocityLagrangianVectorFieldSource.C

\*---------------------------------------------------------------------------*/

#ifndef fanVelocityLagrangianVectorFieldSource_H
#define fanVelocityLagrangianVectorFieldSource_H

#include "LagrangianFieldSources.H"
#include "Function1LagrangianFieldSource.H"
#include "fanDirectionLagrangianVectorFieldSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class fanVelocityLagrangianVectorFieldSource Declaration
\*---------------------------------------------------------------------------*/

class fanVelocityLagrangianVectorFieldSource
:
    public LagrangianVectorFieldSource,
    private Function1LagrangianFieldSource,
    private fanDirectionLagrangianVectorFieldSource
{
private:

    // Private Data

        //- Centreline velocity
        const autoPtr<Function1<vector>> Ucentre_;


    // Private Member Functions

        //- Use the base class' value method
        using Function1LagrangianFieldSource::value;


public:

    //- Runtime type information
    TypeName("fanVelocity");


    // Constructors

        //- Construct from internal field and dictionary
        fanVelocityLagrangianVectorFieldSource
        (
            const regIOobject&,
            const dictionary& dict
        );

        //- Copy constructor setting the internal field reference
        fanVelocityLagrangianVectorFieldSource
        (
            const fanVelocityLagrangianVectorFieldSource&,
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
                new fanVelocityLagrangianVectorFieldSource(*this, iIo)
            );
        }


    //- Destructor
    virtual ~fanVelocityLagrangianVectorFieldSource();


    // Member Functions

        //- Use the base class' value method
        using LagrangianVectorFieldSource::value;

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
