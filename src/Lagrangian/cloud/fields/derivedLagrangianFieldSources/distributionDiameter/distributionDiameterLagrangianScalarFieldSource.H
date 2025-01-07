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
    Foam::distributionDiameterLagrangianScalarFieldSource

Description
    This source condition provides values of diameter randomly sampled from a
    given distribution.

    It is important to use this condition for diameter, rather than the generic
    distribution source condition, as this condition will account for the
    injection of parcels with differing numbers of particles when sampling the
    distribution. If the generic distribution source condition is used when
    injecting parcels of equal mass, for example, then the parcels containing
    many small particles will be over-represented in the sampling. This
    condition corrects for this and produces the correct distribution for
    whatever uniform size the number condition creates.

    This condition must be used with parcel clouds and in conjunction with a
    condition for the number field that specifies a uniform size, such as
    flowRateNumber or totalNumber.

Usage
    \table
        Property      | Description                       | Required? | Default
        distribution  | The distribution                  | yes       |
    \endtable

    Example specification:
    \verbatim
    <LagrangianModelName>
    {
        type            distributionDiameter;
        distribution
        {
            type            RosinRammler;
            Q               0;
            min             0.1 [mm];
            max             1.5 [mm];
            d               1.4 [mm];
            n               2 [];
        }
    }
    \endverbatim

See also
    Foam::distributions
    Foam::flowRateNumberLagrangianScalarFieldSource
    Foam::totalNumberLagrangianScalarFieldSource

SourceFiles
    distributionDiameterLagrangianScalarFieldSource.C

\*---------------------------------------------------------------------------*/

#ifndef distributionDiameterLagrangianScalarFieldSource_H
#define distributionDiameterLagrangianScalarFieldSource_H

#include "LagrangianFieldSources.H"
#include "distribution.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
       Class distributionDiameterLagrangianScalarFieldSource Declaration
\*---------------------------------------------------------------------------*/

class distributionDiameterLagrangianScalarFieldSource
:
    public LagrangianScalarFieldSource
{
private:

    // Private Data

        //- The distribution
        mutable autoPtr<distribution> distribution_;

        //- The time index
        mutable label timeIndex_;

        //- Name of the number field
        const word numberName_;


public:

    //- Runtime type information
    TypeName("distributionDiameter");


    // Constructors

        //- Construct from internal field and dictionary
        distributionDiameterLagrangianScalarFieldSource
        (
            const regIOobject&,
            const dictionary& dict
        );

        //- Copy constructor setting the internal field reference
        distributionDiameterLagrangianScalarFieldSource
        (
            const distributionDiameterLagrangianScalarFieldSource&,
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
                new distributionDiameterLagrangianScalarFieldSource(*this, iIo)
            );
        }


    //- Destructor
    virtual ~distributionDiameterLagrangianScalarFieldSource();


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
