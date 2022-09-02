/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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

#include "thermalBaffleFvPatchScalarField.H"
#include "mappedWallPolyPatch.H"
#include "symmetryPolyPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * /

bool thermalBaffleFvPatchScalarField::primary() const
{
    return patch().boundaryMesh().mesh().name() == polyMesh::defaultRegion;
}


bool thermalBaffleFvPatchScalarField::owner() const
{
    return
        primary()
     && patch().index() < patch().boundaryMesh()[nbrPatch_].index();
}


void thermalBaffleFvPatchScalarField::checkPatches() const
{
    if (!primary()) return;

    const polyPatch& pp = patch().patch();
    const polyPatch& nbrPp = patch().patch().boundaryMesh()[nbrPatch_];

    // The patches should be of mapped type
    auto checkPatchIsMapped = [&](const polyPatch& pp)
    {
        if (!isA<mappedPatchBase>(pp))
        {
            FatalErrorInFunction
                << "Patch field of type \"" << typeName
                << "\" specified for patch \"" << pp.name() << "\" of field \""
                << internalField().name() << "\", but this patch is not of "
                << "type \"" << mappedPatchBase::typeName << "\""
                << exit(FatalError);
        }
    };
    checkPatchIsMapped(pp);
    checkPatchIsMapped(nbrPp);

    const mappedPatchBase& mpp = refCast<const mappedPatchBase>(pp);
    const mappedPatchBase& nbrMpp = refCast<const mappedPatchBase>(nbrPp);

    // The patches should sample a different region
    auto checkPatchMapsDifferentRegion = [&](const mappedPatchBase& mpp)
    {
        if (mpp.sameRegion())
        {
            FatalErrorInFunction
                << "Patch field of type \"" << typeName
                << "\" specified for patch \"" << pp.name() << "\" of field \""
                << internalField().name() << "\", but this patch maps to "
                << "another patch in the same region. It should map to a "
                << "different region; i.e., that of the thermal baffle model."
                << exit(FatalError);
        }
    };
    checkPatchMapsDifferentRegion(mpp);
    checkPatchMapsDifferentRegion(nbrMpp);

    // The sample region of this patch and it's neighbour should be the same,
    // i.e., that of the thermal baffle model
    if (mpp.sampleRegion() != nbrMpp.sampleRegion())
    {
        FatalErrorInFunction
            << "Patch fields of type \"" << typeName
            << "\" specified for patches \"" << pp.name() << "\" and \""
            << nbrPp.name() << "\" of field \"" << internalField().name()
            << "\", but these patches map to different regions \""
            << mpp.sampleRegion() << "\" and \"" << nbrMpp.sampleRegion()
            << ". They should map to the same region; i.e., that of the "
            << "thermal baffle model."
            << exit(FatalError);
    }

    // The sample patch of this patch and it's neighbour should be different,
    // i.e., they should sample opposite ends of the thermal baffle mesh
    if (mpp.samplePatch() == nbrMpp.samplePatch())
    {
        FatalErrorInFunction
            << "Patch fields of type \"" << typeName
            << "\" specified for patches \"" << pp.name() << "\" and \""
            << nbrPp.name() << "\" of field \"" << internalField().name()
            << "\", but these patches map to the same patch; \""
            << mpp.samplePatch() << "\" of region \"" << mpp.sampleRegion()
            << ". They should map to different patches, as these will become "
            << "the patches at opposite ends of the extruded baffle mesh."
            << exit(FatalError);
    }
}


void thermalBaffleFvPatchScalarField::checkPatchFields() const
{
    if (!primary()) return;

    const fvPatch& fvp = patch();
    const fvPatch& nbrFvp = patch().boundaryMesh()[nbrPatch_];

    const fvPatchScalarField& nbrTp =
        nbrFvp.lookupPatchField<volScalarField, scalar>(internalField().name());

    // The neighbour patch field should be of the same type
    if (!isA<thermalBaffleFvPatchScalarField>(nbrTp))
    {
        FatalErrorInFunction
            << "Patch field of type \"" << typeName
            << "\" specified for patch \"" << fvp.name() << "\" of field \""
            << internalField().name() << "\" but the field on the "
            << "neighbouring patch \"" << nbrFvp.name()
            << "\" is of a different type. Both should be of type \""
            << typeName << "\"."
            << exit(FatalError);
    }

    // The neighbour patch field's neighbour should be this patch
    const thermalBaffleFvPatchScalarField& nbrTBp =
        refCast<const thermalBaffleFvPatchScalarField>(nbrTp);
    if (nbrTBp.nbrPatch_ != patch().name())
    {
        FatalErrorInFunction
            << "Patch field of type \"" << typeName
            << "\" on patch \"" << fvp.name() << "\" of field \""
            << internalField().name() << "\" is specified to neighbour "
            << "patch \"" << nbrPatch_ << "\", but this patch does not "
            << "reciprocally neighbour patch \"" << fvp.name() << "\""
            << exit(FatalError);
    }
}


autoPtr<extrudePatchMesh>
thermalBaffleFvPatchScalarField::initBaffleMesh() const
{
    if (!owner())
    {
        FatalErrorInFunction
            << "Baffle mesh is only available to the owner patch in the "
            << "primary region" << exit(FatalError);
    }

    checkPatches();

    const fvMesh& mesh = patch().boundaryMesh().mesh();

    const mappedPatchBase& mpp =
        refCast<const mappedPatchBase>(patch().patch());

    const mappedPatchBase nbrMpp =
        refCast<const mappedPatchBase>
        (patch().patch().boundaryMesh()[nbrPatch_]);

    const List<word> patchNames
    ({
        mpp.samplePatch(),
        nbrMpp.samplePatch(),
        "sides"
    });

    const List<word> patchTypes
    ({
        mappedWallPolyPatch::typeName,
        mappedWallPolyPatch::typeName,
        symmetryPolyPatch::typeName
    });

    List<dictionary> patchDicts(3);
    forAll(patchDicts, patchi)
    {
        patchDicts[patchi].set("nFaces", 0);
        patchDicts[patchi].set("startFace", 0);
    }
    patchDicts[0].add("sampleMode", mpp.sampleModeNames_[mpp.mode()]);
    patchDicts[0].add("sampleRegion", mesh.name());
    patchDicts[0].add("samplePatch", patch().name());
    patchDicts[1].add("sampleMode", mpp.sampleModeNames_[nbrMpp.mode()]);
    patchDicts[1].add("sampleRegion", mesh.name());
    patchDicts[1].add("samplePatch", nbrPatch_);

    List<polyPatch*> patchPtrs(3);
    forAll(patchPtrs, patchi)
    {
        patchPtrs[patchi] = polyPatch::New
        (
            patchTypes[patchi],
            patchNames[patchi],
            patchDicts[patchi],
            patchi,
            mesh.boundaryMesh()
        ).ptr();
    }

    dictionary dict(dict_);
    dict.add("columnCells", false);

    return
        autoPtr<extrudePatchMesh>
        (
            new extrudePatchMesh
            (
                mesh,
                patch(),
                dict,
                mpp.sampleRegion(),
                patchPtrs
            )
        );
}


autoPtr<regionModels::thermalBaffle>
thermalBaffleFvPatchScalarField::initBaffle() const
{
    if (!owner())
    {
        FatalErrorInFunction
            << "Baffle model is only available to the owner patch in the "
            << "primary region" << exit(FatalError);
    }

    checkPatches();

    const fvMesh& mesh = patch().boundaryMesh().mesh();

    const mappedPatchBase& mpp =
        refCast<const mappedPatchBase>(patch().patch());

    dictionary dict(dict_);
    dict.add("regionName", mpp.sampleRegion());

    return autoPtr<regionModels::thermalBaffle>
    (
        new regionModels::thermalBaffle("thermalBaffle", mesh, dict)
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

thermalBaffleFvPatchScalarField::thermalBaffleFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    turbulentTemperatureRadCoupledMixedFvPatchScalarField(p, iF),
    dict_(dictionary::null),
    nbrPatch_(word::null),
    baffleMeshPtr_(),
    bafflePtr_()
{}


thermalBaffleFvPatchScalarField::thermalBaffleFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    turbulentTemperatureRadCoupledMixedFvPatchScalarField(p, iF, dict),
    dict_(dict),
    nbrPatch_(primary() ? dict.lookup<word>("neighbourPatch") : word::null),
    baffleMeshPtr_(owner() ? initBaffleMesh().ptr() : nullptr),
    bafflePtr_(owner() ? initBaffle().ptr() : nullptr)
{}


thermalBaffleFvPatchScalarField::thermalBaffleFvPatchScalarField
(
    const thermalBaffleFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    turbulentTemperatureRadCoupledMixedFvPatchScalarField
    (
        ptf,
        p,
        iF,
        mapper
    ),
    dict_(ptf.dict_),
    baffleMeshPtr_(),
    bafflePtr_()
{}


thermalBaffleFvPatchScalarField::thermalBaffleFvPatchScalarField
(
    const thermalBaffleFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    turbulentTemperatureRadCoupledMixedFvPatchScalarField(ptf, iF),
    dict_(ptf.dict_),
    baffleMeshPtr_(),
    bafflePtr_()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void thermalBaffleFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    turbulentTemperatureRadCoupledMixedFvPatchScalarField::autoMap(m);
}


void thermalBaffleFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    turbulentTemperatureRadCoupledMixedFvPatchScalarField::rmap(ptf, addr);
}


void thermalBaffleFvPatchScalarField::reset
(
    const fvPatchScalarField& ptf
)
{
    turbulentTemperatureRadCoupledMixedFvPatchScalarField::reset(ptf);
}


void thermalBaffleFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    checkPatchFields();

    if (owner())
    {
        bafflePtr_->evolve();
    }

    turbulentTemperatureRadCoupledMixedFvPatchScalarField::updateCoeffs();
}


void thermalBaffleFvPatchScalarField::write(Ostream& os) const
{
    turbulentTemperatureRadCoupledMixedFvPatchScalarField::write(os);

    if (owner())
    {
        forAllConstIter(dictionary, dict_, iter)
        {
            os << *iter;
        }
    }
    else if (primary())
    {
        writeEntry(os, "neighbourPatch", nbrPatch_);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    thermalBaffleFvPatchScalarField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam

// ************************************************************************* //
