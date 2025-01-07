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
    Foam::clouds::shaped

Description
    Base class for clouds with particles with mass

SourceFiles
    shaped.C

\*---------------------------------------------------------------------------*/

#ifndef shaped_H
#define shaped_H

#include "grouped.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace clouds
{

/*---------------------------------------------------------------------------*\
                           Class shaped Declaration
\*---------------------------------------------------------------------------*/

class shaped
{
    // Private Member Functions

        //- Calculate the particles' volumes
        virtual tmp<LagrangianSubScalarField> calcv
        (
            const LagrangianModelRef& model,
            const LagrangianSubMesh& subMesh
        ) const = 0;

        //- Calculate the particles' surface areas
        virtual tmp<LagrangianSubScalarField> calca
        (
            const LagrangianModelRef& model,
            const LagrangianSubMesh& subMesh
        ) const = 0;


public:

    // Public Static Data

        //- Run-time type information
        TypeName("shaped");

        //- Name of volume fields
        static const word vName;


    // Public Data

        //- Volume
        const CloudDerivedField<scalar>& v;

        //- Surface area
        const CloudDerivedField<scalar>& a;

        //- Volume fraction
        const CloudAverageField<scalar>& alpha;


    // Constructors

        //- Construct from a reference to the cloud
        shaped(const cloud&);

        //- Construct from a reference to the cloud and its base classes
        shaped(const cloud&, const grouped&);


    //- Destructor
    virtual ~shaped();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace clouds
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
