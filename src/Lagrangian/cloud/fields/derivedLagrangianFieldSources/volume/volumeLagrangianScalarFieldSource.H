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
    Foam::volumeLagrangianScalarFieldSource

Description
    This source condition provides the volume of introduced particles. It
    converts between the corresponding diameter/etc..., conditions to
    automatically produce a value for the volume. It is selected
    automatically by the functions that need it, so the user should never need
    to specify this condition explicitly.

SourceFiles
    volumeLagrangianScalarFieldSource.C

\*---------------------------------------------------------------------------*/

#ifndef volumeLagrangianScalarFieldSource_H
#define volumeLagrangianScalarFieldSource_H

#include "LagrangianFieldSources.H"
#include "CloudLagrangianFieldSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class volumeLagrangianScalarFieldSource Declaration
\*---------------------------------------------------------------------------*/

class volumeLagrangianScalarFieldSource
:
    public LagrangianScalarFieldSource,
    private CloudLagrangianFieldSource<scalar>
{
public:

    //- Runtime type information
    TypeName("volume");


    // Constructors

        //- Construct from internal field
        volumeLagrangianScalarFieldSource(const regIOobject&);

        //- Construct from internal field and dictionary
        volumeLagrangianScalarFieldSource
        (
            const regIOobject&,
            const dictionary&
        );

        //- Copy constructor setting the internal field reference
        volumeLagrangianScalarFieldSource
        (
            const volumeLagrangianScalarFieldSource&,
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
                new volumeLagrangianScalarFieldSource(*this, iIo)
            );
        }


    //- Destructor
    virtual ~volumeLagrangianScalarFieldSource();


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
