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

\*---------------------------------------------------------------------------*/

#include "gravity.H"
#include "addToRunTimeSelectionTable.H"
#include "lookupUniformDimensionedField.H"
#include "coupledToIncompressibleFluid.H"
#include "coupledToFluid.H"
#include "massive.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace Lagrangian
{
    defineTypeNameAndDebug(gravity, 0);
    addToRunTimeSelectionTable(LagrangianModel, gravity, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Lagrangian::gravity::gravity
(
    const word& name,
    const LagrangianMesh& mesh,
    const dictionary& modelDict,
    const dictionary& stateDict
)
:
    LagrangianModel(name, mesh),
    cloudLagrangianModel(static_cast<const LagrangianModel&>(*this)),
    g(lookupUniformDimensionedField<vector>(mesh.time(), "g"))
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::wordList Foam::Lagrangian::gravity::addSupFields() const
{
    return wordList(1, cloud().U.name());
}


void Foam::Lagrangian::gravity::addSup
(
    const LagrangianSubScalarField& deltaT,
    const LagrangianSubVectorSubField& U,
    LagrangianEqn<vector>& eqn
) const
{
    if
    (
        isCloud<clouds::coupled>()
     && eqn.isPsi(cloud<clouds::coupled>().Uc(U.mesh()))
    )
    {
        // Gravity does not contribute to the carrier momentum equation
    }
    else if (isCloud<clouds::coupledToIncompressibleFluid>())
    {
        const clouds::coupledToIncompressibleFluid& ctifCloud =
            cloud<clouds::coupledToIncompressibleFluid>();

        const dimensionedScalar& rhoByRhoc = ctifCloud.rhoByRhoc;

        eqn.Su += g*(1 - 1/rhoByRhoc);
    }
    else
    {
        eqn.Su += g;
    }
}


void Foam::Lagrangian::gravity::addSup
(
    const LagrangianSubScalarField& deltaT,
    const LagrangianSubScalarSubField& vOrM,
    const LagrangianSubVectorSubField& U,
    LagrangianEqn<vector>& eqn
) const
{
    if
    (
        isCloud<clouds::coupled>()
     && eqn.isPsi(cloud<clouds::coupled>().Uc(U.mesh()))
    )
    {
        // Gravity does not contribute to the carrier momentum equation
    }
    else if (isCloud<clouds::coupledToIncompressibleFluid>())
    {
        const clouds::coupledToIncompressibleFluid& ctifCloud =
            cloud<clouds::coupledToIncompressibleFluid>();

        const dimensionedScalar& rhoByRhoc = ctifCloud.rhoByRhoc;

        eqn.Su += vOrM*g*(1 - 1/rhoByRhoc);
    }
    else if (isCloud<clouds::coupledToFluid>() && isCloud<clouds::massive>())
    {
        const clouds::massive& mCloud = cloud<clouds::massive>();
        const clouds::coupledToFluid& ctfCloud =
            cloud<clouds::coupledToFluid>();

        const LagrangianSubScalarSubField rho(mCloud.rho(U.mesh()));
        const LagrangianSubScalarField& rhoc = ctfCloud.rhoc(U.mesh());

        eqn.Su += vOrM*g*(1 - rhoc/rho);
    }
    else if (isCloud<clouds::coupledToFluid>() && !isCloud<clouds::massive>())
    {
        NotImplemented;
    }
    else
    {
        eqn.Su += vOrM*g;
    }
}


// ************************************************************************* //
