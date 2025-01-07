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
    Foam::flowRateNumberLagrangianScalarFieldSource

Description
    This source condition sets the values of the number field to recover a
    specified volumetric or mass flow-rate. It should be used in conjunction
    with a continuous injection model.

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
        Property           | Description              | Required? | Default
        uniformSize        | The size kept uniform \
                             across all parcels       | yes       |
        volumetricFlowRate | The volumetric flow rate | if massFlowRate \
                                                        is not specified |
        massFlowRate       | The mass flow rate       | if volumetricFlowRate \
                                                        is not specified |
    \endtable

    Example specification:
    \verbatim
    <LagrangianModelName>
    {
        type            flowRateNumber;
        uniformSize     volume; // number, surfaceArea, mass
        volumetricFlowRate 1 [g/s];
    }
    \endverbatim

See also
    Foam::distributionDiameterLagrangianScalarFieldSource

SourceFiles
    flowRateNumberLagrangianScalarFieldSource.C

\*---------------------------------------------------------------------------*/

#ifndef flowRateNumberLagrangianScalarFieldSource_H
#define flowRateNumberLagrangianScalarFieldSource_H

#include "uniformSizeNumberLagrangianScalarFieldSource.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
       Class flowRateNumberLagrangianScalarFieldSource Declaration
\*---------------------------------------------------------------------------*/

class flowRateNumberLagrangianScalarFieldSource
:
    public uniformSizeNumberLagrangianScalarFieldSource
{
private:

    // Private Data

        //- The volumetric flow rate
        autoPtr<Function1<scalar>> volumetricFlowRate_;

        //- The mass flow rate
        autoPtr<Function1<scalar>> massFlowRate_;


public:

    //- Runtime type information
    TypeName("flowRateNumber");


    // Constructors

        //- Construct from internal field and dictionary
        flowRateNumberLagrangianScalarFieldSource
        (
            const regIOobject&,
            const dictionary& dict
        );

        //- Copy constructor setting the internal field reference
        flowRateNumberLagrangianScalarFieldSource
        (
            const flowRateNumberLagrangianScalarFieldSource&,
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
                new flowRateNumberLagrangianScalarFieldSource(*this, iIo)
            );
        }


    //- Destructor
    virtual ~flowRateNumberLagrangianScalarFieldSource();


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
