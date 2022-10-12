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
#include "mappedExtrudedWallPolyPatch.H"
#include "symmetryPolyPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * /

bool Foam::thermalBaffleFvPatchScalarField::primary() const
{
    return patch().boundaryMesh().mesh().name() == polyMesh::defaultRegion;
}


bool Foam::thermalBaffleFvPatchScalarField::owner() const
{
    return
        primary()
     && patch().index() < patch().boundaryMesh()[nbrPatch_].index();
}


void Foam::thermalBaffleFvPatchScalarField::checkPatches() const
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

    const mappedPatchBase& mpb = refCast<const mappedPatchBase>(pp);
    const mappedPatchBase& nbrMpb = refCast<const mappedPatchBase>(nbrPp);

    // The patches should neighbour a different region
    auto checkPatchMapsDifferentRegion = [&](const mappedPatchBase& mpb)
    {
        if (mpb.sameRegion())
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
    checkPatchMapsDifferentRegion(mpb);
    checkPatchMapsDifferentRegion(nbrMpb);

    // The neighbour region of this patch and it's neighbour should be the
    // same, i.e., that of the thermal baffle model
    if (mpb.nbrRegionName() != nbrMpb.nbrRegionName())
    {
        FatalErrorInFunction
            << "Patch fields of type \"" << typeName
            << "\" specified for patches \"" << pp.name() << "\" and \""
            << nbrPp.name() << "\" of field \"" << internalField().name()
            << "\", but these patches map to different regions \""
            << mpb.nbrRegionName() << "\" and \""
            << nbrMpb.nbrRegionName() << ". They should map to the same "
            << "region; i.e., that of the thermal baffle model."
            << exit(FatalError);
    }

    // The neighbour patch of this patch and it's neighbour should be
    // different, i.e., they should map opposite ends of the thermal baffle
    // mesh
    if (mpb.nbrPatchName() == nbrMpb.nbrPatchName())
    {
        FatalErrorInFunction
            << "Patch fields of type \"" << typeName
            << "\" specified for patches \"" << pp.name() << "\" and \""
            << nbrPp.name() << "\" of field \"" << internalField().name()
            << "\", but these patches map to the same patch; \""
            << mpb.nbrPatchName() << "\" of region \""
            << mpb.nbrRegionName() << ". They should map to different "
            << "patches, as these will become the patches at opposite ends of "
            << "the extruded baffle mesh." << exit(FatalError);
    }
}


void Foam::thermalBaffleFvPatchScalarField::checkPatchFields() const
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


Foam::autoPtr<Foam::extrudePatchMesh>
Foam::thermalBaffleFvPatchScalarField::initBaffleMesh() const
{
    if (!owner())
    {
        FatalErrorInFunction
            << "Baffle mesh is only available to the owner patch in the "
            << "primary region" << exit(FatalError);
    }

    checkPatches();

    const fvMesh& mesh = patch().boundaryMesh().mesh();

    const mappedPatchBase& mpb =
        refCast<const mappedPatchBase>(patch().patch());

    const mappedPatchBase nbrMpb =
        refCast<const mappedPatchBase>
        (patch().patch().boundaryMesh()[nbrPatch_]);

    const List<word> patchNames
    ({
        mpb.nbrPatchName(),
        nbrMpb.nbrPatchName(),
        "sides"
    });

    const List<word> patchTypes
    ({
        mappedWallPolyPatch::typeName,
        mappedExtrudedWallPolyPatch::typeName,
        symmetryPolyPatch::typeName
    });

    List<dictionary> patchDicts(3);
    forAll(patchDicts, patchi)
    {
        patchDicts[patchi].set("nFaces", 0);
        patchDicts[patchi].set("startFace", 0);
    }
    patchDicts[0].add("neighbourRegion", mesh.name());
    patchDicts[0].add("neighbourPatch", patch().name());
    patchDicts[1].add("neighbourRegion", mesh.name());
    patchDicts[1].add("neighbourPatch", nbrPatch_);
    patchDicts[1].add("bottomPatch", patchNames[0]);

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
                mpb.nbrRegionName(),
                patchPtrs
            )
        );
}


Foam::autoPtr<Foam::regionModels::thermalBaffle>
Foam::thermalBaffleFvPatchScalarField::initBaffle() const
{
    if (!owner())
    {
        FatalErrorInFunction
            << "Baffle model is only available to the owner patch in the "
            << "primary region" << exit(FatalError);
    }

    checkPatches();

    const fvMesh& mesh = patch().boundaryMesh().mesh();

    const mappedPatchBase& mpb =
        refCast<const mappedPatchBase>(patch().patch());

    dictionary dict(dict_);
    dict.add("regionName", mpb.nbrRegionName());

    return autoPtr<regionModels::thermalBaffle>
    (
        new regionModels::thermalBaffle("thermalBaffle", mesh, dict)
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::thermalBaffleFvPatchScalarField::thermalBaffleFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    coupledTemperatureFvPatchScalarField(p, iF),
    dict_(dictionary::null),
    nbrPatch_(word::null),
    baffleMeshPtr_(),
    bafflePtr_()
{}


Foam::thermalBaffleFvPatchScalarField::thermalBaffleFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    coupledTemperatureFvPatchScalarField(p, iF, dict),
    dict_(dict),
    nbrPatch_(primary() ? dict.lookup<word>("neighbourPatch") : word::null),
    baffleMeshPtr_(owner() ? initBaffleMesh().ptr() : nullptr),
    bafflePtr_(owner() ? initBaffle().ptr() : nullptr)
{}


Foam::thermalBaffleFvPatchScalarField::thermalBaffleFvPatchScalarField
(
    const thermalBaffleFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    coupledTemperatureFvPatchScalarField
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


Foam::thermalBaffleFvPatchScalarField::thermalBaffleFvPatchScalarField
(
    const thermalBaffleFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    coupledTemperatureFvPatchScalarField(ptf, iF),
    dict_(ptf.dict_),
    baffleMeshPtr_(),
    bafflePtr_()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::thermalBaffleFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    coupledTemperatureFvPatchScalarField::autoMap(m);
}


void Foam::thermalBaffleFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    coupledTemperatureFvPatchScalarField::rmap(ptf, addr);
}


void Foam::thermalBaffleFvPatchScalarField::reset
(
    const fvPatchScalarField& ptf
)
{
    coupledTemperatureFvPatchScalarField::reset(ptf);
}


void Foam::thermalBaffleFvPatchScalarField::updateCoeffs()
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

    coupledTemperatureFvPatchScalarField::updateCoeffs();
}


void Foam::thermalBaffleFvPatchScalarField::write(Ostream& os) const
{
    coupledTemperatureFvPatchScalarField::write(os);

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

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        thermalBaffleFvPatchScalarField
    );
}


// ************************************************************************* //
