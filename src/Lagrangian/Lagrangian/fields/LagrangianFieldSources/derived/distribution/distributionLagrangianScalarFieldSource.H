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
    Foam::distributionLagrangianScalarFieldSource

Description
    This source condition provides values of a property randomly sampled from a
    given distribution.

Usage
    \table
        Property          | Description         | Required? | Default
        distribution      | The distribution    | yes       |
    \endtable

    Example specification:
    \verbatim
    <LagrangianModelName>
    {
        type            distribution;
        distribution
        {
            type            normal;
            Q               0;
            min             200 [kg/m^3];
            max             1800 [kg/m^3];
            mu              1000 [kg/m^3];
            sigma           400 [];
        }
    }

See also
    Foam::distributions

SourceFiles
    distributionLagrangianScalarFieldSource.C

\*---------------------------------------------------------------------------*/

#ifndef distributionLagrangianScalarFieldSource_H
#define distributionLagrangianScalarFieldSource_H

#include "LagrangianFieldSources.H"
#include "distribution.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class distributionLagrangianScalarFieldSource Declaration
\*---------------------------------------------------------------------------*/

class distributionLagrangianScalarFieldSource
:
    public LagrangianScalarFieldSource
{
private:

    // Private Data

        //- The distribution
        const autoPtr<distribution> distribution_;

        //- The time index
        mutable label timeIndex_;


public:

    //- Runtime type information
    TypeName("distribution");


    // Constructors

        //- Construct from internal field and dictionary
        distributionLagrangianScalarFieldSource
        (
            const regIOobject&,
            const dictionary& dict
        );

        //- Copy constructor setting the internal field reference
        distributionLagrangianScalarFieldSource
        (
            const distributionLagrangianScalarFieldSource&,
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
                new distributionLagrangianScalarFieldSource(*this, iIo)
            );
        }


    //- Destructor
    virtual ~distributionLagrangianScalarFieldSource();


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
