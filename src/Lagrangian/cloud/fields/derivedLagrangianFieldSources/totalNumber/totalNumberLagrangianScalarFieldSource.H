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
    Foam::totalNumberLagrangianScalarFieldSource

Description
    This source condition sets the values of the number field to recover a
    specified total volume or mass. It should be used in conjunction with an
    instantaneous injection model.

    This condition requires a uniform size to be specified. All injected
    parcels will have the same value of this size. The size can be \c number,
    \c surfaceArea, \c volume or \c mass. This choice relates to the desired
    discretisation of the Lagrangian material. Broadly, the uniform size chosen
    should be that against which the dominant effect of the Lagrangian model
    scales.

    So, if the particles are primarily evaporating droplets of fuel, then the
    effect of each will be primarily a function of the mass of fuel that
    evaporates into the system. The \c uniformSize should therefore be set to
    \c mass.

    But, if the model is purely dynamic, and the dominant effect is viscous
    drag, then each particle will exert a force in proportion with a reference
    area. In this case \c surfaceArea would be a better choice.

    If a distribution of sizes is injected then the choice of uniform size
    interacts with how the distribution is to be sampled. For that reason, it
    is important to use a diameter source condition such as \c
    distributionDiameter which accounts for this, rather than the general \c
    distribution condition.

Usage
    \table
        Property    | Description                | Required? | Default
        uniformSize | The size kept uniform across all parcels | yes |
        volume      | The total volume to inject | if mass is not specified |
        mass        | The total mass to inject   | if volume is not specified |
    \endtable

    Example specification:
    \verbatim
    <LagrangianModelName>
    {
        type            totalNumber;
        uniformSize     volume; // number, surfaceArea, mass
        mass            1 [g];
    }
    \endverbatim

SourceFiles
    totalNumberLagrangianScalarFieldSource.C

\*---------------------------------------------------------------------------*/

#ifndef totalNumberLagrangianScalarFieldSource_H
#define totalNumberLagrangianScalarFieldSource_H

#include "uniformSizeNumberLagrangianScalarFieldSource.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
       Class totalNumberLagrangianScalarFieldSource Declaration
\*---------------------------------------------------------------------------*/

class totalNumberLagrangianScalarFieldSource
:
    public uniformSizeNumberLagrangianScalarFieldSource
{
private:

    // Private Data

        //- Do we have the volume (true) or the mass (false)
        const bool haveVolume_;

        //- The volume or mass
        const scalar volumeOrMass_;


public:

    //- Runtime type information
    TypeName("totalNumber");


    // Constructors

        //- Construct from internal field and dictionary
        totalNumberLagrangianScalarFieldSource
        (
            const regIOobject&,
            const dictionary& dict
        );

        //- Copy constructor setting the internal field reference
        totalNumberLagrangianScalarFieldSource
        (
            const totalNumberLagrangianScalarFieldSource&,
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
                new totalNumberLagrangianScalarFieldSource(*this, iIo)
            );
        }


    //- Destructor
    virtual ~totalNumberLagrangianScalarFieldSource();


    // Member Functions

        //- Return the value for an instantaneous injection
        virtual tmp<LagrangianSubScalarField> value
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
