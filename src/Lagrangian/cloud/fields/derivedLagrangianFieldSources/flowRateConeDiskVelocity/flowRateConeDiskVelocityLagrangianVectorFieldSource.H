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
    Foam::flowRateConeDiskVelocityLagrangianVectorFieldSource

Description
    This source condition provides a conical velocity profile over a disk,
    characterised by inner and outer cone angles. The velocity magnitude is
    obtained from the area of the disk and the flow rate of the corresponding
    flowRateNumber condition. This condition must be used in conjunction with a
    disk injection model and a flow rate condition from the number.

Usage
    \table
        Property      | Description                       | Required? | Default
        thetaInner    | The inner cone angle              | yes       |
        thetaOuter    | The outer cone angle              | yes       |
    \endtable

    Example specification:
    \verbatim
    <LagrangianModelName>
    {
        type            flowRateConeDiskVelocity;
        thetaInner      5 [deg];
        thetaOuter      30 [deg];
    }
    \endverbatim

SourceFiles
    flowRateConeDiskVelocityLagrangianVectorFieldSource.C

\*---------------------------------------------------------------------------*/

#ifndef flowRateConeDiskVelocityLagrangianVectorFieldSource_H
#define flowRateConeDiskVelocityLagrangianVectorFieldSource_H

#include "LagrangianFieldSources.H"
#include "cloudLagrangianFieldSource.H"
#include "coneDiskDirectionLagrangianVectorFieldSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
     Class flowRateConeDiskVelocityLagrangianVectorFieldSource Declaration
\*---------------------------------------------------------------------------*/

class flowRateConeDiskVelocityLagrangianVectorFieldSource
:
    public LagrangianVectorFieldSource,
    private cloudLagrangianFieldSource,
    private coneDiskDirectionLagrangianVectorFieldSource
{
public:

    //- Runtime type information
    TypeName("flowRateConeDiskVelocity");


    // Constructors

        //- Construct from internal field and dictionary
        flowRateConeDiskVelocityLagrangianVectorFieldSource
        (
            const regIOobject&,
            const dictionary& dict
        );

        //- Copy constructor setting the internal field reference
        flowRateConeDiskVelocityLagrangianVectorFieldSource
        (
            const flowRateConeDiskVelocityLagrangianVectorFieldSource&,
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
                new flowRateConeDiskVelocityLagrangianVectorFieldSource
                (
                    *this,
                    iIo
                )
            );
        }


    //- Destructor
    virtual ~flowRateConeDiskVelocityLagrangianVectorFieldSource();


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
