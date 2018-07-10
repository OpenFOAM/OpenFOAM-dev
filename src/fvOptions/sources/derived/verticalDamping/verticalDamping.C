/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2018 OpenFOAM Foundation
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

#include "verticalDamping.H"
#include "fvMesh.H"
#include "fvMatrix.H"
#include "geometricOneField.H"
#include "meshTools.H"
#include "Function1.H"
#include "uniformDimensionedFields.H"
#include "zeroGradientFvPatchField.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(verticalDamping, 0);
    addToRunTimeSelectionTable(option, verticalDamping, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::verticalDamping::add
(
    const volVectorField& alphaRhoU,
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    const uniformDimensionedVectorField& g =
        mesh_.lookupObject<uniformDimensionedVectorField>("g");

    const dimensionedSymmTensor lgg(lambda_*sqr(g)/magSqr(g));

    const DimensionedField<scalar, volMesh>& V = mesh_.V();

    // Calculate the scale
    scalarField s(mesh_.nCells(), ramp_.valid() ? 0 : 1);
    forAll(origins_, i)
    {
        const vectorField& c = mesh_.cellCentres();
        const scalarField x((c - origins_[i]) & directions_[i]);
        s = max(s, ramp_->value(x));
    }

    // Check dimensions
    eqn.dimensions()
      - V.dimensions()*(lgg.dimensions() & alphaRhoU.dimensions());

    // Calculate the force and apply it to the equation
    vectorField force(cells_.size());
    forAll(cells_, i)
    {
        const label c = cells_[i];
        force[i] = V[c]*s[c]*(lgg.value() & alphaRhoU[c]);
    }
    meshTools::constrainDirection(mesh_, mesh_.solutionD(), force);
    forAll(cells_, i)
    {
        const label c = cells_[i];
        eqn.source()[c] += force[i];
    }

    // Write out the force coefficient for debugging
    if (debug && mesh_.time().writeTime())
    {
        volScalarField forceCoeff
        (
            IOobject
            (
                type() + ":forceCoeff",
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            lambda_*mag(g),
            zeroGradientFvPatchField<scalar>::typeName
        );
        forceCoeff.primitiveFieldRef() *= s;
        forceCoeff.write();
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::verticalDamping::verticalDamping
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    cellSetOption(name, modelType, dict, mesh),
    lambda_("lambda", dimless/dimTime, coeffs_.lookup("lambda")),
    ramp_(),
    origins_(),
    directions_()
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::verticalDamping::addSup
(
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    add(eqn.psi(), eqn, fieldi);
}


void Foam::fv::verticalDamping::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    add(rho*eqn.psi(), eqn, fieldi);
}


void Foam::fv::verticalDamping::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    add(alpha*rho*eqn.psi(), eqn, fieldi);
}


bool Foam::fv::verticalDamping::read(const dictionary& dict)
{
    if (cellSetOption::read(dict))
    {
        lambda_ =
            dimensionedScalar
            (
                lambda_.name(),
                lambda_.dimensions(),
                coeffs_.lookup(lambda_.name())
            );

        const bool foundRamp = coeffs_.found("ramp");
        const bool foundOgn = coeffs_.found("origin");
        const bool foundDir = coeffs_.found("direction");
        const bool foundOgns = coeffs_.found("origins");
        const bool foundDirs = coeffs_.found("directions");
        const bool foundAll =
            foundRamp
         && (
                (foundOgn && foundDir && !foundOgns && !foundDirs)
             || (!foundOgn && !foundDir && foundOgns && foundDirs)
            );
         const bool foundAny =
            foundRamp || foundOgn || foundDir || foundOgns || foundDirs;

        if (!foundAll)
        {
            ramp_ = autoPtr<Function1<scalar>>();
            origins_.clear();
            directions_.clear();
        }

        if (foundAll)
        {
            ramp_ = Function1<scalar>::New("ramp", coeffs_);
            if (foundOgn)
            {
                origins_.setSize(1);
                directions_.setSize(1);
                coeffs_.lookup("origin") >> origins_.last();
                coeffs_.lookup("direction") >> directions_.last();
            }
            else
            {
                coeffs_.lookup("origins") >> origins_;
                coeffs_.lookup("directions") >> directions_;

                if
                (
                    origins_.size() == 0
                 || directions_.size() == 0
                 || origins_.size() != directions_.size()
                )
                {
                    FatalErrorInFunction
                        << "The same, non-zero number of origins and "
                        << "directions must be provided" << exit(FatalError);
                }
            }
            forAll(directions_, i)
            {
                directions_[i] /= mag(directions_[i]);
            }
        }

        if (!foundAll && foundAny)
        {
            WarningInFunction
                << "The ramping specification is incomplete. \"ramp\", "
                << "\"origin\" and \"direction\" (or \"origins\" and "
                << "\"directions\"), must all be specified in order to ramp "
                << "the damping. The damping will be applied uniformly across "
                << "the cell set." << endl << endl;
        }

        fieldNames_ = wordList(1, coeffs_.lookupOrDefault<word>("U", "U"));

        applied_.setSize(1, false);

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
