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
    Foam::clouds::coupledToFluid

Description
    Base class for clouds which are coupled to a fluid

SourceFiles
    coupledToFluid.C

\*---------------------------------------------------------------------------*/

#ifndef coupledToFluid_H
#define coupledToFluid_H

#include "coupled.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace clouds
{

/*---------------------------------------------------------------------------*\
                           Class coupledToFluid Declaration
\*---------------------------------------------------------------------------*/

class coupledToFluid
:
    public coupled
{
private:

    // Private Member Data

        //- Reference to the mesh
        const LagrangianMesh& mesh_;

        //- Cached carrier density field
        tmp<volScalarField> trhocVf_;


    // Private Member Functions

        //- Get the carrier density field
        virtual tmp<volScalarField> getRhocVf() const;

        //- Calculate the carrier kinematic viscosity
        virtual tmp<LagrangianSubScalarField> calcNuc
        (
            const LagrangianModelRef& model,
            const LagrangianSubMesh& subMesh
        ) const;


protected:

    // Protected Member Functions

        //- Update the carrier density
        void updateRhoc();


public:

    // Public Static Data

        //- Run-time type information
        TypeName("coupledToFluid");


    // Public Data

        //- Carrier density
        const CarrierField<scalar>& rhoc;

        //- Carrier dynamic viscosity
        const CarrierField<scalar>& muc;


    // Constructors

        //- Construct from a reference to the cloud
        coupledToFluid(const cloud&);


    //- Destructor
    virtual ~coupledToFluid();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace clouds
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
