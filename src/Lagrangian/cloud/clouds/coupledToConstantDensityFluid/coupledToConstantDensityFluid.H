/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025-2026 OpenFOAM Foundation
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
    Foam::clouds::coupledToConstantDensityFluid

Description
    Base class for clouds which are coupled to a constant density fluid

SourceFiles
    coupledToConstantDensityFluid.C

\*---------------------------------------------------------------------------*/

#ifndef coupledToConstantDensityFluid_H
#define coupledToConstantDensityFluid_H

#include "coupled.H"
#include "physicalProperties.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace clouds
{

/*---------------------------------------------------------------------------*\
                 Class coupledToConstantDensityFluid Declaration
\*---------------------------------------------------------------------------*/

class coupledToConstantDensityFluid
:
    public coupled
{
private:

    // Private Member Data

        //- Reference to the mesh
        const LagrangianMesh& mesh_;

        //- Physical properties dictionary
        mutable autoPtr<Foam::physicalProperties> physicalPropertiesPtr_;

        //- Reference to the constant cloud density. Only valid if the
        //  Eulerian is multiphase. Null if the Eulerian is single phase.
        const dimensionedScalar& rho_;

        //- Reference to the constant carrier density. Only valid if the
        //  Eulerian is multiphase. Null if the Eulerian is single phase.
        const dimensionedScalar& rhoc_;


    // Private Member Functions

        //- Get the physical properties dictionary
        const Foam::physicalProperties& physicalProperties() const;

        //- Return a uniform field of values equal to one. Used in volume
        //  transfers to represent the Eulerian volume per unit volume.
        static tmp<LagrangianSubScalarField> one
        (
            const LagrangianSubMesh&
        );


protected:

    // Protected Member Functions

        //- A uniform field of values equal to one. Used in volume
        //  transfers to represent the Eulerian volume per unit volume.
        const CloudDerivedField<scalar>& onec;

        //- A uniform field of values equal to one. Used in volume
        //  transfers to represent the Eulerian volume per unit volume.
        const CloudDerivedField<scalar>& onecPhase;


public:

    // Public Static Data

        //- Run-time type information
        TypeName("coupledToConstantDensityFluid");


    // Public Data

        //- Cloud/carrier density ratio
        const dimensionedScalar rhoByRhoc;


    // Constructors

        //- Construct from a reference to the cloud
        coupledToConstantDensityFluid(const cloud& c, const dictionary& dict);


    //- Destructor
    virtual ~coupledToConstantDensityFluid();


    // Member Functions

        //- Cloud density
        const dimensionedScalar& rho() const;

        //- Carrier density
        const dimensionedScalar& rhoc() const;

        //- Corresponding Eulerian phase density
        const dimensionedScalar& rhocPhase() const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace clouds
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
