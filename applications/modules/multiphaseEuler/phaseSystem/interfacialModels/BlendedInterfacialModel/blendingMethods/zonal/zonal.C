/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2026 OpenFOAM Foundation
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

#include "zonal.H"
#include "zonalBoundaryIndices.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace blendingMethods
{
    defineTypeNameAndDebug(zonal, 0);
    addToRunTimeSelectionTable(blendingMethod, zonal, dictionary);
}
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::blendingMethods::zonal::f
(
    const tmp<volScalarField>& defaultF,
    const PtrList<volScalarField>& zoneFs
) const
{
    const fvMesh& mesh = defaultF().mesh();

    const boundaryIndices& zoneIndexBf =
        boundaryIndices::New(mesh, zoneIndices_);

    tmp<volScalarField> tf(defaultF);
    volScalarField& f = tf.ref();
    volScalarField::BoundaryField& fBf = f.boundaryFieldRef();

    defaultF.clear();

    forAll(zoneIndices_, zonei)
    {
        const cellZone& cells = mesh.cellZones()[zoneIndices_[zonei]];

        UIndirectList<scalar>(f.primitiveFieldRef(), cells) =
            UIndirectList<scalar>(zoneFs[zonei].primitiveField(), cells);
    }

    forAll(mesh.boundary(), patchi)
    {
        forAll(mesh.boundary()[patchi], patchFacei)
        {
            const label zonei = zoneIndexBf[patchi][patchFacei];

            if (zonei != -1)
            {
                fBf[patchi][patchFacei] =
                    zoneFs[zonei].boundaryField()[patchi][patchFacei];
            }
        }
    }

    return tf;
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * /

Foam::tmp<Foam::volScalarField> Foam::blendingMethods::zonal::fContinuous
(
    const UPtrList<const volScalarField>& alphas,
    const label index
) const
{
    tmp<volScalarField> defaultFContinuous =
        blendingMethod::fContinuous
        (
            defaultBlending_(),
            alphas,
            index
        );

    PtrList<volScalarField> zoneFContinuouss(zoneIndices_.size());
    forAll(zoneIndices_, zonei)
    {
        zoneFContinuouss.set
        (
            zonei,
            blendingMethod::fContinuous
            (
                zoneBlendings_[zonei],
                alphas,
                index
            ).ptr()
        );
    }

    return f(defaultFContinuous, zoneFContinuouss);
}


Foam::tmp<Foam::volScalarField> Foam::blendingMethods::zonal::fDisplaced
(
    const UPtrList<const volScalarField>& alphas,
    const label displacingPhasei
) const
{
    tmp<volScalarField> defaultFDisplaced =
        blendingMethod::fDisplaced
        (
            defaultBlending_(),
            alphas,
            displacingPhasei
        );

    PtrList<volScalarField> zoneFDisplaceds(zoneIndices_.size());
    forAll(zoneIndices_, zonei)
    {
        zoneFDisplaceds.set
        (
            zonei,
            blendingMethod::fDisplaced
            (
                zoneBlendings_[zonei],
                alphas,
                displacingPhasei
            ).ptr()
        );
    }

    return f(defaultFDisplaced, zoneFDisplaceds);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::blendingMethods::zonal::zonal
(
    const dictionary& dict,
    const phaseInterface& interface,
    const bool allowDisplaced
)
:
    blendingMethod(interface),
    defaultBlending_
    (
        blendingMethod::New
        (
            "default",
            dict.subDict("default"),
            interface,
            allowDisplaced
        )
    ),
    zoneIndices_(),
    zoneBlendings_()
{
    const dictionary& blendingsDict = dict.subDict("zones");

    zoneIndices_.resize(blendingsDict.size());
    zoneBlendings_.resize(blendingsDict.size());

    label zonei = 0;

    forAllConstIter(dictionary, blendingsDict, iter)
    {
        zoneIndices_[zonei] =
            interface.mesh().cellZones().findIndex(iter().keyword());

        zoneBlendings_.set
        (
            zonei,
            blendingMethod::New
            (
                iter().keyword(),
                iter().dict(),
                interface
            )
        );

        zonei ++;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::blendingMethods::zonal::~zonal()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::blendingMethods::zonal::canBeContinuous(const label index) const
{
    if (defaultBlending_->canBeContinuous(index))
    {
        return true;
    }

    forAll(zoneIndices_, zonei)
    {
        if (zoneBlendings_[zonei].canBeContinuous(index))
        {
            return true;
        }
    }

    return false;
}


bool Foam::blendingMethods::zonal::canSegregate() const
{
    if (defaultBlending_->canSegregate())
    {
        return true;
    }

    forAll(zoneIndices_, zonei)
    {
        if (zoneBlendings_[zonei].canSegregate())
        {
            return true;
        }
    }

    return false;
}


bool Foam::blendingMethods::zonal::isDisplacedBy
(
    const label displacingPhasei
) const
{
    if (defaultBlending_->isDisplacedBy(displacingPhasei))
    {
        return true;
    }

    forAll(zoneIndices_, zonei)
    {
        if (zoneBlendings_[zonei].isDisplacedBy(displacingPhasei))
        {
            return true;
        }
    }

    return false;
}


bool Foam::blendingMethods::zonal::functionOfAlphas() const
{
    return false;
}


// ************************************************************************* //
