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
    Foam::noneStateLagrangianLabelFieldSource

Description
    This source condition provides a label corresponding to the
    LagrangianState::none enumeration.

SourceFiles
    noneStateLagrangianLabelFieldSource.C

\*---------------------------------------------------------------------------*/

#ifndef noneStateLagrangianLabelFieldSource_H
#define noneStateLagrangianLabelFieldSource_H

#include "LagrangianFieldSources.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class noneStateLagrangianLabelFieldSource Declaration
\*---------------------------------------------------------------------------*/

class noneStateLagrangianLabelFieldSource
:
    public LagrangianLabelFieldSource
{
public:

    //- Runtime type information
    TypeName("noneState");


    // Constructors

        //- Inherit constructors
        using LagrangianLabelFieldSource::LagrangianFieldSource;

        //- Construct and return a clone setting the internal field reference
        virtual autoPtr<LagrangianLabelFieldSource> clone
        (
            const regIOobject& iIo
        ) const
        {
            return autoPtr<LagrangianLabelFieldSource>
            (
                new noneStateLagrangianLabelFieldSource(*this, iIo)
            );
        }


    //- Destructor
    virtual ~noneStateLagrangianLabelFieldSource();


    // Member Functions

        //- Return the value for an instantaneous injection
        virtual tmp<LagrangianSubLabelField> value
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
