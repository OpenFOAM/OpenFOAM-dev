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
    Foam::specificHeatCapacityLagrangianScalarFieldSource

Description
    Source condition for specific heat capacity. Applied to the specific heat
    capacity field automatically by the Lagrangian thermodynamic model. This
    should not need to be specified manually.

SourceFiles
    specificHeatCapacityLagrangianScalarFieldSource.C

\*---------------------------------------------------------------------------*/

#ifndef specificHeatCapacityLagrangianScalarFieldSource_H
#define specificHeatCapacityLagrangianScalarFieldSource_H

#include "LagrangianFieldSources.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
Class specificHeatCapacityLagrangianScalarFieldSource Declaration
\*---------------------------------------------------------------------------*/

class specificHeatCapacityLagrangianScalarFieldSource
:
    public LagrangianScalarFieldSource
{
public:

    //- Runtime type information
    TypeName("specificHeatCapacity");


    // Constructors

        //- Inherit constructors
        using LagrangianScalarFieldSource::LagrangianScalarFieldSource;

        //- Construct and return a clone setting the internal field reference
        virtual autoPtr<LagrangianScalarFieldSource> clone
        (
            const regIOobject& iIo
        ) const
        {
            return autoPtr<LagrangianScalarFieldSource>
            (
                new specificHeatCapacityLagrangianScalarFieldSource(*this, iIo)
            );
        }


    //- Destructor
    virtual ~specificHeatCapacityLagrangianScalarFieldSource();


    // Member Functions

        //- Return the value for an instantaneous injection
        virtual tmp<LagrangianSubScalarField> value
        (
            const LagrangianInjection&,
            const LagrangianSubMesh&
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
