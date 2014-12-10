/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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

#include "pressureGradientExplicitSource.H"
#include "fvMatrices.H"
#include "DimensionedField.H"
#include "IFstream.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(pressureGradientExplicitSource, 0);

    addToRunTimeSelectionTable
    (
        option,
        pressureGradientExplicitSource,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::pressureGradientExplicitSource::writeProps
(
    const scalar gradP
) const
{
    // Only write on output time
    if (mesh_.time().outputTime())
    {
        IOdictionary propsDict
        (
            IOobject
            (
                name_ + "Properties",
                mesh_.time().timeName(),
                "uniform",
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            )
        );
        propsDict.add("gradient", gradP);
        propsDict.regIOobject::write();
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::pressureGradientExplicitSource::pressureGradientExplicitSource
(
    const word& sourceName,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    option(sourceName, modelType, dict, mesh),
    Ubar_(coeffs_.lookup("Ubar")),
    gradP0_(0.0),
    dGradP_(0.0),
    flowDir_(Ubar_/mag(Ubar_)),
    invAPtr_(NULL)
{
    coeffs_.lookup("fieldNames") >> fieldNames_;

    if (fieldNames_.size() != 1)
    {
        FatalErrorIn
        (
            "Foam::fv::pressureGradientExplicitSource::"
            "pressureGradientExplicitSource"
            "("
                "const word&, "
                "const word&, "
                "const dictionary&, "
                "const fvMesh&"
            ")"
        )   << "Source can only be applied to a single field.  Current "
            << "settings are:" << fieldNames_ << exit(FatalError);
    }

    applied_.setSize(fieldNames_.size(), false);

    // Read the initial pressure gradient from file if it exists
    IFstream propsFile
    (
        mesh_.time().timePath()/"uniform"/(name_ + "Properties")
    );

    if (propsFile.good())
    {
        Info<< "    Reading pressure gradient from file" << endl;
        dictionary propsDict(dictionary::null, propsFile);
        propsDict.lookup("gradient") >> gradP0_;
    }

    Info<< "    Initial pressure gradient = " << gradP0_ << nl << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::pressureGradientExplicitSource::correct(volVectorField& U)
{
    const scalarField& rAU = invAPtr_().internalField();

    // Integrate flow variables over cell set
    scalar magUbarAve = 0.0;
    scalar rAUave = 0.0;
    const scalarField& cv = mesh_.V();
    forAll(cells_, i)
    {
        label cellI = cells_[i];
        scalar volCell = cv[cellI];
        magUbarAve += (flowDir_ & U[cellI])*volCell;
        rAUave += rAU[cellI]*volCell;
    }

    // Collect across all processors
    reduce(magUbarAve, sumOp<scalar>());
    reduce(rAUave, sumOp<scalar>());

    // Volume averages
    magUbarAve /= V_;
    rAUave /= V_;

    // Calculate the pressure gradient increment needed to adjust the average
    // flow-rate to the desired value
    dGradP_ = (mag(Ubar_) - magUbarAve)/rAUave;

    // Apply correction to velocity field
    forAll(cells_, i)
    {
        label cellI = cells_[i];
        U[cellI] += flowDir_*rAU[cellI]*dGradP_;
    }

    scalar gradP = gradP0_ + dGradP_;

    Info<< "Pressure gradient source: uncorrected Ubar = " << magUbarAve
        << ", pressure gradient = " << gradP << endl;

    writeProps(gradP);
}


void Foam::fv::pressureGradientExplicitSource::addSup
(
    fvMatrix<vector>& eqn,
    const label fieldI
)
{
    DimensionedField<vector, volMesh> Su
    (
        IOobject
        (
            name_ + fieldNames_[fieldI] + "Sup",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector("zero", eqn.dimensions()/dimVolume, vector::zero)
    );

    scalar gradP = gradP0_ + dGradP_;

    UIndirectList<vector>(Su, cells_) = flowDir_*gradP;

    eqn += Su;
}


void Foam::fv::pressureGradientExplicitSource::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldI
)
{
    this->addSup(eqn, fieldI);
}


void Foam::fv::pressureGradientExplicitSource::setValue
(
    fvMatrix<vector>& eqn,
    const label
)
{
    if (invAPtr_.empty())
    {
        invAPtr_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    name_ + ":invA",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                1.0/eqn.A()
            )
        );
    }
    else
    {
        invAPtr_() = 1.0/eqn.A();
    }

    gradP0_ += dGradP_;
    dGradP_ = 0.0;
}


// ************************************************************************* //
