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

#include "rigidBodySectionalForcesBase.H"
#include "rigidBodyMotion.H"
#include "motionSolver_fvMeshMover.H"
#include "motionSolver.H"
#include "transformField.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

vector AxCrossX(const tensor& A, const symmTensor& xSqr)
{
    return
        vector
        (
            A.yx()*xSqr.zx() + A.yy()*xSqr.zy() + A.yz()*xSqr.zz()
          - A.zx()*xSqr.yx() - A.zy()*xSqr.yy() - A.zz()*xSqr.yz(),

            A.zx()*xSqr.xx() + A.zy()*xSqr.xy() + A.zz()*xSqr.xz()
          - A.xx()*xSqr.zx() - A.xy()*xSqr.zy() - A.xz()*xSqr.zz(),

            A.xx()*xSqr.yx() + A.xy()*xSqr.yy() + A.xz()*xSqr.yz()
          - A.yx()*xSqr.xx() - A.yy()*xSqr.xy() - A.yz()*xSqr.xz()
        );
}

tmp<vectorField> AxCrossX(const tensor& A, const symmTensorField& xSqr)
{
    tmp<vectorField> tresult(new vectorField(xSqr.size()));
    vectorField& result = tresult.ref();
    forAll(result, i)
    {
        result[i] = AxCrossX(A, xSqr[i]);
    }
    return tresult;
}

}


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(rigidBodySectionalForcesBase, 0);
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

const Foam::RBD::rigidBodyMotion&
Foam::functionObjects::rigidBodySectionalForcesBase::motion() const
{
    const fvMeshMovers::motionSolver& mover =
        refCast<const fvMeshMovers::motionSolver>(mesh_.mover());

    return refCast<const RBD::rigidBodyMotion>(mover.motion());
}


Foam::label
Foam::functionObjects::rigidBodySectionalForcesBase::bodyIndex() const
{
    return motion().bodyIndex(body_);
}


const Foam::RBD::rigidBody&
Foam::functionObjects::rigidBodySectionalForcesBase::body() const
{
    return motion().bodies()[bodyIndex()];
}


Foam::label Foam::functionObjects::rigidBodySectionalForcesBase::axisi() const
{
    return axisi_;
}


const Foam::word&
Foam::functionObjects::rigidBodySectionalForcesBase::axisName() const
{
    return coordSet::axisTypeNames_[axis_];
}


Foam::vector Foam::functionObjects::rigidBodySectionalForcesBase::axis() const
{
    vector result = vector::zero;
    result[axisi_] = 1;
    return result;
}


const Foam::point&
Foam::functionObjects::rigidBodySectionalForcesBase::localOrigin() const
{
    return localOrigin_;
}


Foam::vector Foam::functionObjects::rigidBodySectionalForcesBase::normal() const
{
    return motion().X0(bodyIndex()).E().vectorComponent(axisi_);
}


Foam::point Foam::functionObjects::rigidBodySectionalForcesBase::origin() const
{
    return motion().X0(bodyIndex()).inv().transformPoint(localOrigin_);
}


void Foam::functionObjects::rigidBodySectionalForcesBase::addFluid
(
    vectorField& force,
    vectorField& moment
) const
{
    vectorField dForce(force.size(), vector::zero);
    vectorField dMoment(moment.size(), vector::zero);

    // Calculate the fluid contribution
    sectionalForcesBase::addFluid(dForce, dMoment);

    // Rotate into the body's coordinate frame
    const tensor E = motion().X0(bodyIndex()).E();
    transform(dForce, E, dForce);
    transform(dMoment, E, dMoment);

    // Add to the total
    force += dForce;
    moment += dMoment;
}


void Foam::functionObjects::rigidBodySectionalForcesBase::addBody
(
    vectorField& force,
    vectorField& moment
) const
{
    tmp<scalarField> tdistances = this->distances();
    const scalarField& distances = tdistances();

    // Get the velocity and acceleration in the frame of the body
    const vector& v0 = motion().v(bodyIndex()).l();
    const vector& w0 = motion().v(bodyIndex()).w();
    const vector vDot0 = motion().a(bodyIndex()).l() - (v0 ^ w0);
    const vector wDot0 = motion().a(bodyIndex()).w();
    const tensor wDotXPlusWXWX(*wDot0 + (*w0 & *w0));

    // Get the moments of the body sections
    const scalarField mu0s
    (
        body().sectionMu0s(axisi_, distances + localOrigin_[axisi_])
    );
    const vectorField mu1s
    (
        body().sectionMu1s(axisi_, distances + localOrigin_[axisi_])
    );
    const symmTensorField mu2s
    (
        body().sectionMu2s(axisi_, distances + localOrigin_[axisi_])
    );

    // Calculate the body forces and moments on the intervals
    const vectorField intervalForces
    (
        - (vDot0*mu0s)
        - (wDotXPlusWXWX & mu1s)
    );
    const vectorField intervalMoments
    (
          (vDot0 ^ mu1s)
        + AxCrossX(wDotXPlusWXWX, mu2s)
        - (localOrigin_ ^ intervalForces)
    );

    // Cumulatively sum the interval forces to obtain the sectional forces
    vector f = vector::zero, m = vector::zero;
    forAllReverse(intervalForces, i)
    {
        f += intervalForces[i];
        m += intervalMoments[i];
        force[i] += f;
        moment[i] += m - (distances[i]*axis() ^ f);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::rigidBodySectionalForcesBase::
rigidBodySectionalForcesBase
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    sectionalForcesBase(name, runTime, dict),
    body_()
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::rigidBodySectionalForcesBase::
~rigidBodySectionalForcesBase()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::rigidBodySectionalForcesBase::read
(
    const dictionary& dict
)
{
    sectionalForcesBase::read(dict);

    body_ = dict.lookup<word>("body");

    localOrigin_ = dict.lookupOrDefault<point>("origin", body().c());

    axis_ = coordSet::axisTypeNames_.read(dict.lookup("axis"));

    switch (coordSet::axisTypeNames_.read(dict.lookup("axis")))
    {
        case coordSet::axisType::X:
        case coordSet::axisType::Y:
        case coordSet::axisType::Z:
            break;
        default:
            FatalIOErrorInFunction(dict)
                << "Axis for " << typeName << " function '" << name()
                << " must be "
                << coordSet::axisTypeNames_[coordSet::axisType::X] << ", "
                << coordSet::axisTypeNames_[coordSet::axisType::Y] << " or "
                << coordSet::axisTypeNames_[coordSet::axisType::Z]
                << exit(FatalIOError);
    }

    axisi_ =
        static_cast<label>(axis_)
      - static_cast<label>(coordSet::axisType::X);

    return true;
}


void Foam::functionObjects::rigidBodySectionalForcesBase::movePoints
(
    const polyMesh& mesh
)
{
    if (&mesh == &mesh_)
    {
        // Don't clear the weights and directions. The patch transformation
        // should be solid-body, so these don't change. The actual patch
        // geometry does, though, so re-generate that.
        clearPatchGeom();
    }
}


// ************************************************************************* //
