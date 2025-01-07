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
    Foam::clouds::spherical

Description
    Base class for clouds with spherical particles

SourceFiles
    spherical.C

\*---------------------------------------------------------------------------*/

#ifndef spherical_H
#define spherical_H

#include "shaped.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace clouds
{

/*---------------------------------------------------------------------------*\
                    Class spherical Declaration
\*---------------------------------------------------------------------------*/

class spherical
:
    public shaped
{
private:

    // Private Member Functions

        //- Return the IO object for the diameter field
        static IOobject dIo(const cloud&);

        //- Calculate the particles' volumes
        virtual tmp<LagrangianSubScalarField> calcv
        (
            const LagrangianModelRef& model,
            const LagrangianSubMesh& subMesh
        ) const;

        //- Calculate the particles' surface areas
        virtual tmp<LagrangianSubScalarField> calca
        (
            const LagrangianModelRef& model,
            const LagrangianSubMesh& subMesh
        ) const;


protected:

    // Protected Member Functions

        //- Correct the shape to match the given volume
        void correct(const LagrangianSubScalarSubField& v);


public:

    // Public Static Data

        //- Run-time type information
        TypeName("spherical");


    // Public Data

        //- Diameter
        CloudStateField<scalar>& d;


    // Constructors

        //- Construct from a reference to the cloud
        spherical(const cloud&);

        //- Construct from a reference to the cloud and its base classes
        spherical(const cloud&, const grouped&);


    //- Destructor
    virtual ~spherical();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace clouds
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
