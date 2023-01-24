/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2023 OpenFOAM Foundation
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

#include "alphatPhaseChangeWallFunctionBase.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "phaseInterface.H"
#include "phaseSystem.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(alphatPhaseChangeWallFunctionBase, 0);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

alphatPhaseChangeWallFunctionBase::alphatPhaseChangeWallFunctionBase()
:
    phaseName_(word::null),
    otherPhaseName_(word::null)
{}


alphatPhaseChangeWallFunctionBase::alphatPhaseChangeWallFunctionBase
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    phaseName_(iF.group()),
    otherPhaseName_(dict.lookup("otherPhase"))
{
    if (phaseName_== otherPhaseName_)
    {
        FatalErrorInFunction
            << "otherPhase should be the name of the vapour phase that "
            << "corresponds to the liquid base or vice versa" << nl
            << "This phase: " << phaseName_ << nl
            << "otherPhase: " << otherPhaseName_
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

alphatPhaseChangeWallFunctionBase::~alphatPhaseChangeWallFunctionBase()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool alphatPhaseChangeWallFunctionBase::activeInterface
(
    const phaseInterface& interface
) const
{
    const phaseSystem& fluid = interface.fluid();

    return
        interface.contains(fluid.phases()[phaseName_])
     && interface.contains(fluid.phases()[otherPhaseName_]);
}


void alphatPhaseChangeWallFunctionBase::write(Ostream& os) const
{
    writeEntry(os, "otherPhase", otherPhaseName_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam

// ************************************************************************* //
