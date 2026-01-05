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
    Foam::compressibilityLagrangianScalarFieldSource

Description
    Source condition for compressibility. Applied to the compressibility field
    automatically by the Lagrangian thermodynamic model. This should not need
    to be specified manually.

SourceFiles
    compressibilityLagrangianScalarFieldSource.C

\*---------------------------------------------------------------------------*/

#ifndef compressibilityLagrangianScalarFieldSource_H
#define compressibilityLagrangianScalarFieldSource_H

#include "LagrangianFieldSources.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
Class compressibilityLagrangianScalarFieldSource Declaration
\*---------------------------------------------------------------------------*/

class compressibilityLagrangianScalarFieldSource
:
    public LagrangianScalarFieldSource
{
public:

    //- Runtime type information
    TypeName("compressibility");


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
                new compressibilityLagrangianScalarFieldSource(*this, iIo)
            );
        }


    //- Destructor
    virtual ~compressibilityLagrangianScalarFieldSource();


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
