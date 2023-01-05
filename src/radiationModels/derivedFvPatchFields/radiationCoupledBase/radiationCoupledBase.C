/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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

#include "radiationCoupledBase.H"
#include "volFields.H"
#include "mappedPatchBase.H"
#include "fvPatchFieldMapper.H"
#include "radiationModel.H"
#include "opaqueSolid.H"
#include "absorptionEmissionModel.H"

// * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(radiationCoupledBase, 0);
}

template<>
const char* Foam::NamedEnum
<
    Foam::radiationCoupledBase::emissivityMethodType,
    2
>::names[] =
{
    "solidRadiation",
    "lookup"
};

const Foam::NamedEnum
<
    Foam::radiationCoupledBase::emissivityMethodType,
    2
> Foam::radiationCoupledBase::emissivityMethodTypeNames_;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiationCoupledBase::radiationCoupledBase
(
    const fvPatch& patch,
    const word& calculationType,
    const scalarField& emissivity
)
:
    patch_(patch),
    method_(emissivityMethodTypeNames_[calculationType]),
    emissivity_(emissivity)
{}


Foam::radiationCoupledBase::radiationCoupledBase
(
    const fvPatch& patch,
    const word& calculationType,
    const scalarField& emissivity,
    const fvPatchFieldMapper& mapper
)
:
    patch_(patch),
    method_(emissivityMethodTypeNames_[calculationType]),
    emissivity_(mapper(emissivity))
{}


Foam::radiationCoupledBase::radiationCoupledBase
(
    const fvPatch& patch,
    const dictionary& dict
)
:
    patch_(patch),
    method_(emissivityMethodTypeNames_.read(dict.lookup("emissivityMode")))
{
    switch (method_)
    {
        case SOLIDRADIATION:
        {
            emissivity_ = scalarField(patch_.size(), scalar(0));
            break;
        }

        case LOOKUP:
        {
            emissivity_ = scalarField("emissivity", dict, patch_.size());
            break;
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor    * * * * * * * * * * * * * * //

Foam::radiationCoupledBase::~radiationCoupledBase()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField> Foam::radiationCoupledBase::emissivity() const
{
    switch (method_)
    {
        case SOLIDRADIATION:
        {
            // Get the coupling information from the mappedPatchBase
            const mappedPatchBase& mpp =
                mappedPatchBase::getMap(patch_.patch());

            const polyMesh& nbrMesh = mpp.nbrMesh();

            const radiationModels::opaqueSolid& radiation =
                nbrMesh.lookupObject<radiationModels::opaqueSolid>
                (
                    "radiationProperties"
                );

            const fvMesh& nbrFvMesh = refCast<const fvMesh>(nbrMesh);

            const fvPatch& nbrPatch =
                nbrFvMesh.boundary()[mpp.nbrPolyPatch().index()];

            // NOTE: for an opaqueSolid the absorptionEmission model returns the
            // emissivity of the surface rather than the emission coefficient
            // and the input specification MUST correspond to this.
            return
                mpp.distribute
                (
                    radiation.absorptionEmission().e()().boundaryField()
                    [
                        nbrPatch.index()
                    ]
                );
        }
        break;

        case LOOKUP:
        {
            // Return local value
            return emissivity_;
        }

        default:
        {
            FatalErrorInFunction
                << "Unimplemented method " << method_ << endl
                << "Please set 'emissivity' to one of "
                << emissivityMethodTypeNames_.toc()
                << exit(FatalError);
        }
        break;
    }

    return scalarField(0);
}


void Foam::radiationCoupledBase::autoMap(const fvPatchFieldMapper& m)
{
    m(emissivity_, emissivity_);
}


void Foam::radiationCoupledBase::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    const radiationCoupledBase& mrptf =
        refCast<const radiationCoupledBase>(ptf);

    emissivity_.rmap(mrptf.emissivity_, addr);
}


void Foam::radiationCoupledBase::reset
(
    const fvPatchScalarField& ptf
)
{
    const radiationCoupledBase& mrptf =
        refCast<const radiationCoupledBase>(ptf);

    emissivity_.reset(mrptf.emissivity_);
}


void Foam::radiationCoupledBase::write(Ostream& os) const
{
    writeEntry(os, "emissivityMode", emissivityMethodTypeNames_[method_]);
    writeEntry(os, "emissivity", emissivity_);
}


// ************************************************************************* //
