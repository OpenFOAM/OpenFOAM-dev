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

#include "bXiKernelCorr.H"
#include "bXiIgnition.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
    namespace fv
    {
        defineTypeNameAndDebug(bXiKernelCorr, 0);

        addToRunTimeSelectionTable
        (
            fvModel,
            bXiKernelCorr,
            dictionary
        );
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::bXiKernelCorr::readCoeffs(const dictionary& dict)
{
    duration_.read(dict, mesh().time().userUnits());
    bMin_.readIfPresent(dict);
}


template<class Type>
inline auto Foam::fv::bXiKernelCorr::bFunc(const Type& b, const Type& R) const
{
    // Analytical solution of b-Xi model with linear Xi and heat release
    return (R - (R - 1)*b)*(1 - b);

    // Fit to Gaussian kernel position distribution
    // return pow(max(1 - b, 0.0), 0.8)/pow(max(b, bMin_.value()), 0.2);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::bXiKernelCorr::bXiKernelCorr
(
    const word& name,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    fvModel(name, modelType, mesh, dict),
    zone_(mesh, coeffs(dict)),
    kernelShape_(kernelShape::New(mesh, coeffs(dict))),
    ignition_
    (
        refCast<const fv::bXiIgnition>
        (
            fvModels::New(mesh)[coeffs(dict).lookup<word>("ignition")]
        )
    ),
    duration_("duration", mesh.time().userUnits(), coeffs(dict)),
    bMin_
    (
        "bMin",
        dimless,
        coeffs(dict),
        mesh.lookupObject<solvers::XiFluid>(solver::typeName).bMin()/2
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::wordList Foam::fv::bXiKernelCorr::addSupFields() const
{
    return ignition_.addSupFields();
}


void Foam::fv::bXiKernelCorr::addSup
(
    const volScalarField& rho,
    const volScalarField& b,
    fvMatrix<scalar>& eqn
) const
{
    if (!ignition_.igniting(duration_)) return;

    const labelList& cells = zone_.zone();

    const scalarField bc(b, cells);
    const scalar bMin = returnReduce(min(bc), minOp<scalar>());

    // Check if min(b) is < 1, i.e. ignition has occurred
    // and min(b) > small, i.e. burnout has not occurred requiring correction
    if (bMin < 1 - bMin_.value() && bMin > bMin_.value())
    {
        if (debug)
        {
            Info<< type()
                << ": applying source to " << eqn.psi().name() << endl;
        }

        const volScalarField& rhou = mesh().lookupObject<volScalarField>
        (
            IOobject::groupName("rho", "u")
        );

        const volScalarField& rhob = mesh().lookupObject<volScalarField>
        (
            IOobject::groupName("rho", "b")
        );

        const volScalarField& mgb = mesh().lookupObject<volScalarField>("mgb");
        const volScalarField& Sl = mesh().lookupObject<volScalarField>("Su");
        const volScalarField& Xi = mesh().lookupObject<volScalarField>("Xi");

        const scalarField& V = mesh().V();

        const scalarField Vc(V, cells);
        const scalarField Rc(scalarField(rhou, cells)/scalarField(rhob, cells));

        // Calculate volume of kernel
        const dimensionedScalar Vk("Vk", dimVolume, gSum((1 - bc)*Vc));

        // Calculate kernel area from its volume
        const scalar Ak(kernelShape_->Ak(Vk).value());

        // Calculate volume integral of kernel distribution function
        const scalar Vkd = gSum(bFunc(bc, Rc)*bc*Vc);

        scalarField& Sp = eqn.diag();

        forAll(cells, i)
        {
            const label celli = cells[i];
            const scalar b = bc[i];

            // Add kernel propagation correction source
            Sp[celli] -=
                Vc[i]*rhou[celli]*Sl[celli]*Xi[celli]
               *max(Ak*bFunc(b, Rc[i])/Vkd - mgb[celli]/max(b, bMin), 0);
        }
    }
}


void Foam::fv::bXiKernelCorr::topoChange
(
    const polyTopoChangeMap& map
)
{
    zone_.topoChange(map);
}


void Foam::fv::bXiKernelCorr::mapMesh(const polyMeshMap& map)
{
    zone_.mapMesh(map);
}


void Foam::fv::bXiKernelCorr::distribute
(
    const polyDistributionMap& map
)
{
    zone_.distribute(map);
}


bool Foam::fv::bXiKernelCorr::movePoints()
{
    zone_.movePoints();
    return true;
}


bool Foam::fv::bXiKernelCorr::read(const dictionary& dict)
{
    if (fvModel::read(dict))
    {
        zone_.read(coeffs(dict));
        readCoeffs(coeffs(dict));
        return true;
    }
    else
    {
        return false;
    }

    return false;
}


// ************************************************************************* //
