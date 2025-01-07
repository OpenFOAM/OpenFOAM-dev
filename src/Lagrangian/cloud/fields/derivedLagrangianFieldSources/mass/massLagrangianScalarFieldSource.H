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
    Foam::massLagrangianScalarFieldSource

Description
    This source condition provides the mass of introduced particles. It
    converts between the corresponding volume/diameter/density/etc...,
    conditions to automatically produce a value for the mass. It is selected
    automatically by the clouds that need it, so the user should never need
    to specify this condition explicitly.

SourceFiles
    massLagrangianScalarFieldSource.C

\*---------------------------------------------------------------------------*/

#ifndef massLagrangianScalarFieldSource_H
#define massLagrangianScalarFieldSource_H

#include "LagrangianFieldSources.H"
#include "CloudLagrangianFieldSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class massLagrangianScalarFieldSource Declaration
\*---------------------------------------------------------------------------*/

class massLagrangianScalarFieldSource
:
    public LagrangianScalarFieldSource,
    private CloudLagrangianFieldSource<scalar>
{
public:

    //- Runtime type information
    TypeName("mass");


    // Constructors

        //- Construct from internal field
        massLagrangianScalarFieldSource(const regIOobject&);

        //- Construct from internal field and dictionary
        massLagrangianScalarFieldSource(const regIOobject&, const dictionary&);

        //- Copy constructor setting the internal field reference
        massLagrangianScalarFieldSource
        (
            const massLagrangianScalarFieldSource&,
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
                new massLagrangianScalarFieldSource(*this, iIo)
            );
        }


    //- Destructor
    virtual ~massLagrangianScalarFieldSource();


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
