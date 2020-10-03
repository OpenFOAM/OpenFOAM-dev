/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2020 OpenFOAM Foundation
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

#include "alphatPhaseChangeWallFunctionFvPatchScalarField.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(alphatPhaseChangeWallFunctionFvPatchScalarField, 0);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

alphatPhaseChangeWallFunctionFvPatchScalarField::
alphatPhaseChangeWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    alphatPhaseJayatillekeWallFunctionFvPatchScalarField(p, iF),
    otherPhaseName_(word::null),
    relax_(1),
    dmdtf_(p.size(), 0)
{}


alphatPhaseChangeWallFunctionFvPatchScalarField::
alphatPhaseChangeWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    alphatPhaseJayatillekeWallFunctionFvPatchScalarField(p, iF, dict),
    otherPhaseName_(dict.lookup("otherPhase")),
    relax_(dict.lookupOrDefault<scalar>("relax", 1)),
    dmdtf_(p.size(), 0)
{
    // Check that otherPhaseName != this phase
    if (internalField().group() == otherPhaseName_)
    {
        FatalErrorInFunction
            << "otherPhase should be the name of the vapor phase that "
            << "corresponds to the liquid base or vice versa" << nl
            << "This phase: " << internalField().group() << nl
            << "otherPhase: " << otherPhaseName_
            << abort(FatalError);
    }

    if (dict.found("dmdtf"))
    {
        dmdtf_ = scalarField("dmdtf", dict, p.size());
    }
}


alphatPhaseChangeWallFunctionFvPatchScalarField::
alphatPhaseChangeWallFunctionFvPatchScalarField
(
    const alphatPhaseChangeWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    alphatPhaseJayatillekeWallFunctionFvPatchScalarField(ptf, p, iF, mapper),
    otherPhaseName_(ptf.otherPhaseName_),
    relax_(ptf.relax_),
    dmdtf_(mapper(ptf.dmdtf_))
{}


alphatPhaseChangeWallFunctionFvPatchScalarField::
alphatPhaseChangeWallFunctionFvPatchScalarField
(
    const alphatPhaseChangeWallFunctionFvPatchScalarField& awfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    alphatPhaseJayatillekeWallFunctionFvPatchScalarField(awfpsf, iF),
    otherPhaseName_(awfpsf.otherPhaseName_),
    relax_(awfpsf.relax_),
    dmdtf_(awfpsf.dmdtf_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool alphatPhaseChangeWallFunctionFvPatchScalarField::
activePhasePair(const phasePairKey& phasePair) const
{
    if (phasePair == phasePairKey(otherPhaseName_, internalField().group()))
    {
        return true;
    }
    else
    {
        return false;
    }
}


const scalarField&
alphatPhaseChangeWallFunctionFvPatchScalarField::dmdtf() const
{
    return dmdtf_;
}


const scalarField& alphatPhaseChangeWallFunctionFvPatchScalarField::
dmdtf(const phasePairKey& phasePair) const
{
    if (activePhasePair(phasePair))
    {
        return dmdtf_;
    }
    else
    {
        FatalErrorInFunction
            << " dmdtf requested for invalid phasePair!"
            << abort(FatalError);

        return dmdtf_;
    }
}


void alphatPhaseChangeWallFunctionFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    alphatPhaseJayatillekeWallFunctionFvPatchScalarField::autoMap(m);

    m(dmdtf_, dmdtf_);
}


void alphatPhaseChangeWallFunctionFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    alphatPhaseJayatillekeWallFunctionFvPatchScalarField::rmap(ptf, addr);

    const alphatPhaseChangeWallFunctionFvPatchScalarField& tiptf =
        refCast<const alphatPhaseChangeWallFunctionFvPatchScalarField>(ptf);

    dmdtf_.rmap(tiptf.dmdtf_, addr);
}


void alphatPhaseChangeWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    alphatPhaseJayatillekeWallFunctionFvPatchScalarField::write(os);

    writeEntry(os, "otherPhase", otherPhaseName_);
    writeEntry(os, "relax", relax_);
    writeEntry(os, "dmdtf", dmdtf_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam

// ************************************************************************* //
