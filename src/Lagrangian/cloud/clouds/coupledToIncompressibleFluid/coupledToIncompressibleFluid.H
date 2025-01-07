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
    Foam::clouds::coupledToIncompressibleFluid

Description
    Base class for clouds which are coupled to an incompressible fluid

SourceFiles
    coupledToIncompressibleFluid.C

\*---------------------------------------------------------------------------*/

#ifndef coupledToIncompressibleFluid_H
#define coupledToIncompressibleFluid_H

#include "coupled.H"
#include "physicalProperties.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace clouds
{

/*---------------------------------------------------------------------------*\
                 Class coupledToIncompressibleFluid Declaration
\*---------------------------------------------------------------------------*/

class coupledToIncompressibleFluid
:
    public coupled
{
private:

    // Private Member Data

        //- Reference to the mesh
        const LagrangianMesh& mesh_;

        //- Physical properties dictionary
        Foam::physicalProperties physicalProperties_;

        //- Cached carrier kinematic viscosity field
        tmp<volScalarField> tnucVf_;


    // Private Member Functions

        //- Get the carrier kinematic viscosity field
        virtual tmp<volScalarField> getNucVf() const;

        //- Calculate the carrier kinematic viscosity
        virtual tmp<LagrangianSubScalarField> calcNuc
        (
            const LagrangianModelRef& model,
            const LagrangianSubMesh& subMesh
        ) const;


protected:

    // Protected Member Functions

        //- Update the carrier kinematic viscosity
        void updateNuc();


public:

    // Public Static Data

        //- Run-time type information
        TypeName("coupledToIncompressibleFluid");


    // Public Data

        //- Cloud/carrier density ratio
        const dimensionedScalar rhoByRhoc;

        //- Carrier kinematic viscosity
        const CarrierField<scalar>& nuc;


    // Constructors

        //- Construct from a reference to the cloud
        coupledToIncompressibleFluid(const cloud&);


    //- Destructor
    virtual ~coupledToIncompressibleFluid();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace clouds
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
