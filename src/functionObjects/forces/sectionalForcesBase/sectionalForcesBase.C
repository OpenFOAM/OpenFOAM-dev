/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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

#include "indirectPrimitivePatch.H"
#include "sectionalForcesBase.H"
#include "surfaceFields.H"
#include "volFields.H"
#include "writeFile.H"

#include "fvcGrad.H"
#include "surfaceInterpolate.H"

#include "incompressibleMomentumTransportModel.H"
#include "compressibleMomentumTransportModel.H"
#include "phaseIncompressibleMomentumTransportModel.H"
#include "phaseCompressibleMomentumTransportModel.H"

#include "polyTopoChangeMap.H"
#include "polyMeshMap.H"
#include "polyDistributionMap.H"

#include "forces.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(sectionalForcesBase, 0);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::SurfaceField<Type>>
Foam::functionObjects::sectionalForcesBase::timesAlpha
(
    const tmp<SurfaceField<Type>> psi
) const
{
    if (phaseName_ != word::null)
    {
        const volScalarField& alpha =
            mesh().lookupObject<volScalarField>
            (
                IOobject::groupName("alpha", phaseName_)
            );

        psi.ref() *= fvc::interpolate(alpha);
    }

    return psi;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::functionObjects::sectionalForcesBase::timesAlpha
(
    const tmp<Field<Type>> psi,
    const label patchi
) const
{
    if (phaseName_ != word::null)
    {
        const volScalarField& alpha =
            mesh().lookupObject<volScalarField>
            (
                IOobject::groupName("alpha", phaseName_)
            );

        psi.ref() *= alpha.boundaryField()[patchi];
    }

    return psi;
}


template<class Type>
Foam::tmp<Foam::SurfaceField<Type>>
Foam::functionObjects::sectionalForcesBase::timesRho
(
    const tmp<SurfaceField<Type>> psi
) const
{
    if (rhoName_ != "rhoInf")
    {
        const volScalarField & rho =
            mesh().lookupObject<volScalarField>(rhoName_);

        psi.ref() *= fvc::interpolate(rho);
    }
    else
    {
        psi.ref() *= dimensionedScalar(dimDensity, rhoRef_);
    }

    return psi;
}


template<class Type>
Foam::tmp<Foam::SurfaceField<Type>>
Foam::functionObjects::sectionalForcesBase::timesAlphaRho
(
    const tmp<SurfaceField<Type>> psi
) const
{
    if (phaseName_ != word::null && rhoName_ != "rhoInf")
    {
        const volScalarField& alpha =
            mesh().lookupObject<volScalarField>
            (
                IOobject::groupName("alpha", phaseName_)
            );

        const volScalarField & rho =
            mesh().lookupObject<volScalarField>(rhoName_);

        psi.ref() *= fvc::interpolate(alpha*rho);

        return psi;
    }
    else
    {
        return timesRho(timesAlpha(psi));
    }
}


Foam::tmp<Foam::volScalarField>
Foam::functionObjects::sectionalForcesBase::p() const
{
    const volScalarField& p = obr_.lookupObject<volScalarField>(pName_);

    if (p.dimensions() == dimPressure)
    {
        return p;
    }
    else if (p.dimensions() == dimPressure/dimDensity)
    {
        if (rhoName_ != "rhoInf")
        {
            FatalErrorInFunction
                << "kinematic pressure found but no 'rhoInf' specified"
                << exit(FatalError);
        }

        return dimensionedScalar(dimDensity, rhoRef_)*p;
    }
    else
    {
        FatalErrorInFunction
            << "pressure dimensions not recognised"
            << exit(FatalError);

        return tmp<volScalarField>(nullptr);
    }
}


Foam::tmp<Foam::surfaceVectorField>
Foam::functionObjects::sectionalForcesBase::devTau() const
{
    typedef incompressible::momentumTransportModel icoModel;
    typedef compressible::momentumTransportModel cmpModel;
    typedef phaseIncompressible::momentumTransportModel phaseIcoModel;
    typedef phaseCompressible::momentumTransportModel phaseCmpModel;

    const word& modelName = momentumTransportModel::typeName;
    const word phaseModelName = IOobject::groupName(modelName, phaseName_);

    if (obr_.foundObject<icoModel>(modelName))
    {
        const incompressible::momentumTransportModel& model =
            obr_.lookupObject<icoModel>(modelName);

        return timesAlphaRho(model.devSigma());
    }
    else if (obr_.foundObject<cmpModel>(modelName))
    {
        const cmpModel& model =
            obr_.lookupObject<cmpModel>(modelName);

        return timesAlpha(model.devTau());
    }
    else if (obr_.foundObject<phaseIcoModel>(phaseModelName))
    {
        const phaseIcoModel& model =
            obr_.lookupObject<phaseIcoModel>(phaseModelName);

        return timesRho(model.devSigma());
    }
    else if (obr_.foundObject<phaseCmpModel>(phaseModelName))
    {
        const phaseCmpModel& model =
            obr_.lookupObject<phaseCmpModel>(phaseModelName);

        return model.devTau();
    }
    else
    {
        FatalErrorInFunction
            << "No valid model for viscous stress calculation"
            << exit(FatalError);

        return surfaceVectorField::null();
    }
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::sectionalForcesBase::clear()
{
    weightsPtr_.clear();
}


const Foam::indirectPrimitivePatch&
Foam::functionObjects::sectionalForcesBase::patch() const
{
    if (!patchPtr_.valid())
    {
        labelList patchFaces;
        forAllConstIter(labelHashSet, patchSet_, iter)
        {
            const label ppi = iter.key();
            const polyPatch& pp = mesh().boundaryMesh()[ppi];

            patchFaces.append(identityMap(pp.size()) + pp.start());
        }

        patchPtr_.reset
        (
            new indirectPrimitivePatch
            (
                IndirectList<face>(mesh().faces(), patchFaces),
                mesh().points()
            )
        );
    }

    return patchPtr_();
}


void Foam::functionObjects::sectionalForcesBase::clearPatch()
{
    if (patchPtr_.valid())
    {
        patchPtr_.clear();
    }
}


void Foam::functionObjects::sectionalForcesBase::clearPatchGeom()
{
    if (patchPtr_.valid())
    {
        patchPtr_->clearGeom();
    }
}


Foam::tmp<Foam::scalarField>
Foam::functionObjects::sectionalForcesBase::patchPointDistances() const
{
    return (patch().localPoints() - origin()) & normal();
}


const Foam::List<Foam::patchCutPlot::weight>&
Foam::functionObjects::sectionalForcesBase::weights() const
{
    if (!weightsPtr_.valid())
    {
        weightsPtr_.reset
        (
            new List<patchCutPlot::weight>
            (
                patchCutPlot::calcWeights
                (
                    patch(),
                    patchPointDistances(),
                    distances(),
                    false,
                    false
                )
            )
        );
    }

    return weightsPtr_();
}


Foam::fileName Foam::functionObjects::sectionalForcesBase::outputPath() const
{
    return
        time_.globalPath()
       /writeFile::outputPrefix
       /(mesh_.name() != polyMesh::defaultRegion ? mesh_.name() : word())
       /name()
       /time_.name();
}


void Foam::functionObjects::sectionalForcesBase::addFluid
(
    vectorField& force,
    vectorField& moment
) const
{
    tmp<scalarField> tdistances = this->distances();
    const scalarField& distances = tdistances();

    // Get the pressure
    tmp<volScalarField> tp = this->p();
    const volScalarField& p = tp();

    // Get the stress tensor
    tmp<surfaceVectorField> tdevTau = this->devTau();
    const surfaceVectorField& devTau = tdevTau();

    // Compute patch-face fluid forces and moments around the origin
    label patchFacei = 0;
    vectorField patchFaceForces(patch().size());
    vectorField patchFaceMoments(patch().size());
    forAllConstIter(labelHashSet, patchSet_, iter)
    {
        const label ppi = iter.key();
        const polyPatch& pp = mesh().boundaryMesh()[ppi];

        const vectorField f
        (
            timesAlpha
            (
                pp.faceNormals()*(p.boundaryField()[ppi] - pRef_)
              + devTau.boundaryField()[ppi],
                ppi
            )
        );

        SubList<vector>(patchFaceForces, pp.size(), patchFacei) =
            f;
        SubList<vector>(patchFaceMoments, pp.size(), patchFacei) =
            (pp.faceCentres() - origin()) ^ f;

        patchFacei += pp.size();
    }

    // Construct the total fluid forces on the intervals
    const vectorField intervalForces
    (
        cutPlot::applyWeights
        (
            distances.size() - 1,
            weights(),
            patchFaceForces
        )
    );
    const vectorField intervalMoments
    (
        cutPlot::applyWeights
        (
            distances.size() - 1,
            weights(),
            patchFaceMoments
        )
    );

    // Cumulatively sum the interval forces to obtain the sectional forces
    vector f = vector::zero, m = vector::zero;
    forAllReverse(intervalForces, i)
    {
        f += intervalForces[i];
        m += intervalMoments[i];
        force[i] += f;
        moment[i] += m - (distances[i]*normal() ^ f);
    }

    // Check consistency with the forces function object if the calculation
    // range spans the entire patch
    if
    (
        debug
     && distances.first() < gMin(patchPointDistances())
     && gMax(patchPointDistances()) < distances.last()
    )
    {
        wordList patchNames;
        forAllConstIter(labelHashSet, patchSet_, iter)
        {
            const label ppi = iter.key();
            const polyPatch& pp = mesh().boundaryMesh()[ppi];

            patchNames.append(pp.name());
        }

        functionObjects::forces forcesFunctionObject
        (
            functionObjects::forces::typeName,
            mesh().time(),
            dictionary::entries
            (
                "type", functionObjects::forces::typeName,
                "patches", patchNames,
                "CofR", origin() + distances[0]*normal()
            )
        );

        forcesFunctionObject.calcForcesMoments();

        Info<< functionObject::typeName << "s::"
            << functionObjects::forces::typeName << ":" << nl
            << "    force = " << forcesFunctionObject.forceEff() << nl
            << "    moment = " << forcesFunctionObject.momentEff() << nl
            << functionObject::typeName << "s::"
            << type() << ":" << nl
            << "    force = " << force[0] << nl
            << "    moment = " << moment[0] << nl << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::sectionalForcesBase::sectionalForcesBase
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    patchSet_(),
    patchPtr_(nullptr),
    pName_(word::null),
    UName_(word::null),
    rhoName_(word::null),
    phaseName_(word::null),
    rhoRef_(NaN),
    pRef_(NaN),
    weightsPtr_(nullptr)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::sectionalForcesBase::~sectionalForcesBase()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::sectionalForcesBase::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    patchSet_ = mesh().boundaryMesh().patchSet(dict);

    // Optional phase entry
    phaseName_ = dict.lookupOrDefault<word>("phase", word::null);

    // Optional p, U, and rho entries
    pName_ =
        dict.lookupOrDefault<word>
        (
            "p",
            IOobject::groupName("p", phaseName_)
        );
    UName_ =
        dict.lookupOrDefault<word>
        (
            "U",
            IOobject::groupName("U", phaseName_)
        );
    rhoName_ =
        dict.lookupOrDefault<word>
        (
            "rho",
            IOobject::groupName("rho", phaseName_)
        );

    // Reference density needed for incompressible calculations
    if (rhoName_ == "rhoInf")
    {
        dict.lookup("rhoInf") >> rhoRef_;
    }

    // Reference pressure, 0 by default
    pRef_ = dict.lookupOrDefault<scalar>("pRef", 0.0);

    return true;
}


Foam::wordList Foam::functionObjects::sectionalForcesBase::fields() const
{
    return wordList::null();
}


bool Foam::functionObjects::sectionalForcesBase::execute()
{
    return true;
}


bool Foam::functionObjects::sectionalForcesBase::end()
{
    return true;
}


void Foam::functionObjects::sectionalForcesBase::movePoints
(
    const polyMesh& mesh
)
{
    if (&mesh == &mesh_)
    {
        clear();
        clearPatchGeom();
    }
}


void Foam::functionObjects::sectionalForcesBase::topoChange
(
    const polyTopoChangeMap& map
)
{
    if (&map.mesh() == &mesh_)
    {
        clear();
        clearPatch();
    }
}


void Foam::functionObjects::sectionalForcesBase::mapMesh
(
    const polyMeshMap& map
)
{
    if (&map.mesh() == &mesh_)
    {
        clear();
        clearPatch();
    }
}


void Foam::functionObjects::sectionalForcesBase::distribute
(
    const polyDistributionMap& map
)
{
    if (&map.mesh() == &mesh_)
    {
        clear();
        clearPatch();
    }
}


// ************************************************************************* //
