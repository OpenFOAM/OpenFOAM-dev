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

#include "rotorDiskSource.H"
#include "addToRunTimeSelectionTable.H"
#include "trimModel.H"
#include "fvMatrices.H"
#include "geometricOneField.H"
#include "syncTools.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
    namespace fv
    {
        defineTypeNameAndDebug(rotorDiskSource, 0);
        addToRunTimeSelectionTable(option, rotorDiskSource, dictionary);
    }

    template<> const char* NamedEnum<fv::rotorDiskSource::geometryModeType, 2>::
        names[] =
    {
        "auto",
        "specified"
    };

    const NamedEnum<fv::rotorDiskSource::geometryModeType, 2>
        fv::rotorDiskSource::geometryModeTypeNames_;

    template<> const char* NamedEnum<fv::rotorDiskSource::inletFlowType, 3>::
        names[] =
    {
        "fixed",
        "surfaceNormal",
        "local"
    };

    const NamedEnum<fv::rotorDiskSource::inletFlowType, 3>
        fv::rotorDiskSource::inletFlowTypeNames_;
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::fv::rotorDiskSource::checkData()
{
    // set inflow type
    switch (selectionMode())
    {
        case smCellSet:
        case smCellZone:
        case smAll:
        {
            // set the profile ID for each blade section
            profiles_.connectBlades(blade_.profileName(), blade_.profileID());
            switch (inletFlow_)
            {
                case ifFixed:
                {
                    coeffs_.lookup("inletVelocity") >> inletVelocity_;
                    break;
                }
                case ifSurfaceNormal:
                {
                    scalar UIn
                    (
                        readScalar(coeffs_.lookup("inletNormalVelocity"))
                    );
                    inletVelocity_ = -coordSys_.R().e3()*UIn;
                    break;
                }
                case ifLocal:
                {
                    // do nothing
                    break;
                }
                default:
                {
                    FatalErrorIn("void rotorDiskSource::checkData()")
                        << "Unknown inlet velocity type" << abort(FatalError);
                }
            }


            break;
        }
        default:
        {
            FatalErrorIn("void rotorDiskSource::checkData()")
                << "Source cannot be used with '"
                << selectionModeTypeNames_[selectionMode()]
                << "' mode.  Please use one of: " << nl
                << selectionModeTypeNames_[smCellSet] << nl
                << selectionModeTypeNames_[smCellZone] << nl
                << selectionModeTypeNames_[smAll]
                << exit(FatalError);
        }
    }
}


void Foam::fv::rotorDiskSource::setFaceArea(vector& axis, const bool correct)
{
    area_ = 0.0;

    static const scalar tol = 0.8;

    const label nInternalFaces = mesh_.nInternalFaces();
    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();
    const vectorField& Sf = mesh_.Sf();
    const scalarField& magSf = mesh_.magSf();

    vector n = vector::zero;

    // calculate cell addressing for selected cells
    labelList cellAddr(mesh_.nCells(), -1);
    UIndirectList<label>(cellAddr, cells_) = identity(cells_.size());
    labelList nbrFaceCellAddr(mesh_.nFaces() - nInternalFaces, -1);
    forAll(pbm, patchI)
    {
        const polyPatch& pp = pbm[patchI];

        if (pp.coupled())
        {
            forAll(pp, i)
            {
                label faceI = pp.start() + i;
                label nbrFaceI = faceI - nInternalFaces;
                label own = mesh_.faceOwner()[faceI];
                nbrFaceCellAddr[nbrFaceI] = cellAddr[own];
            }
        }
    }

    // correct for parallel running
    syncTools::swapBoundaryFaceList(mesh_, nbrFaceCellAddr);

    // add internal field contributions
    for (label faceI = 0; faceI < nInternalFaces; faceI++)
    {
        const label own = cellAddr[mesh_.faceOwner()[faceI]];
        const label nbr = cellAddr[mesh_.faceNeighbour()[faceI]];

        if ((own != -1) && (nbr == -1))
        {
            vector nf = Sf[faceI]/magSf[faceI];

            if ((nf & axis) > tol)
            {
                area_[own] += magSf[faceI];
                n += Sf[faceI];
            }
        }
        else if ((own == -1) && (nbr != -1))
        {
            vector nf = Sf[faceI]/magSf[faceI];

            if ((-nf & axis) > tol)
            {
                area_[nbr] += magSf[faceI];
                n -= Sf[faceI];
            }
        }
    }


    // add boundary contributions
    forAll(pbm, patchI)
    {
        const polyPatch& pp = pbm[patchI];
        const vectorField& Sfp = mesh_.Sf().boundaryField()[patchI];
        const scalarField& magSfp = mesh_.magSf().boundaryField()[patchI];

        if (pp.coupled())
        {
            forAll(pp, j)
            {
                const label faceI = pp.start() + j;
                const label own = cellAddr[mesh_.faceOwner()[faceI]];
                const label nbr = nbrFaceCellAddr[faceI - nInternalFaces];
                const vector nf = Sfp[j]/magSfp[j];

                if ((own != -1) && (nbr == -1) && ((nf & axis) > tol))
                {
                    area_[own] += magSfp[j];
                    n += Sfp[j];
                }
            }
        }
        else
        {
            forAll(pp, j)
            {
                const label faceI = pp.start() + j;
                const label own = cellAddr[mesh_.faceOwner()[faceI]];
                const vector nf = Sfp[j]/magSfp[j];

                if ((own != -1) && ((nf & axis) > tol))
                {
                    area_[own] += magSfp[j];
                    n += Sfp[j];
                }
            }
        }
    }

    if (correct)
    {
        reduce(n, sumOp<vector>());
        axis = n/mag(n);
    }

    if (debug)
    {
        volScalarField area
        (
            IOobject
            (
                name_ + ":area",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("0", dimArea, 0)
        );
        UIndirectList<scalar>(area.internalField(), cells_) = area_;

        Info<< type() << ": " << name_ << " writing field " << area.name()
            << endl;

        area.write();
    }
}


void Foam::fv::rotorDiskSource::createCoordinateSystem()
{
    // construct the local rotor co-prdinate system
    vector origin(vector::zero);
    vector axis(vector::zero);
    vector refDir(vector::zero);

    geometryModeType gm =
        geometryModeTypeNames_.read(coeffs_.lookup("geometryMode"));

    switch (gm)
    {
        case gmAuto:
        {
            // determine rotation origin (cell volume weighted)
            scalar sumV = 0.0;
            const scalarField& V = mesh_.V();
            const vectorField& C = mesh_.C();
            forAll(cells_, i)
            {
                const label cellI = cells_[i];
                sumV += V[cellI];
                origin += V[cellI]*C[cellI];
            }
            reduce(origin, sumOp<vector>());
            reduce(sumV, sumOp<scalar>());
            origin /= sumV;

            // determine first radial vector
            vector dx1(vector::zero);
            scalar magR = -GREAT;
            forAll(cells_, i)
            {
                const label cellI = cells_[i];
                vector test = C[cellI] - origin;
                if (mag(test) > magR)
                {
                    dx1 = test;
                    magR = mag(test);
                }
            }
            reduce(dx1, maxMagSqrOp<vector>());
            magR = mag(dx1);

            // determine second radial vector and cross to determine axis
            forAll(cells_, i)
            {
                const label cellI = cells_[i];
                vector dx2 = C[cellI] - origin;
                if (mag(dx2) > 0.5*magR)
                {
                    axis = dx1 ^ dx2;
                    if (mag(axis) > SMALL)
                    {
                        break;
                    }
                }
            }
            reduce(axis, maxMagSqrOp<vector>());
            axis /= mag(axis);

            // correct the axis direction using a point above the rotor
            {
                vector pointAbove(coeffs_.lookup("pointAbove"));
                vector dir = pointAbove - origin;
                dir /= mag(dir);
                if ((dir & axis) < 0)
                {
                    axis *= -1.0;
                }
            }

            coeffs_.lookup("refDirection") >> refDir;

            localAxesRotation_.reset
            (
                new localAxesRotation
                (
                    mesh_,
                    axis,
                    origin,
                    cells_
                )
            );

            // set the face areas and apply correction to calculated axis
            // e.g. if cellZone is more than a single layer in thickness
            setFaceArea(axis, true);

            break;
        }
        case gmSpecified:
        {
            coeffs_.lookup("origin") >> origin;
            coeffs_.lookup("axis") >> axis;
            coeffs_.lookup("refDirection") >> refDir;

            localAxesRotation_.reset
            (
                new localAxesRotation
                (
                    mesh_,
                    axis,
                    origin,
                    cells_
                )
            );

            setFaceArea(axis, false);

            break;
        }
        default:
        {
            FatalErrorIn("rotorDiskSource::createCoordinateSystem()")
                << "Unknown geometryMode " << geometryModeTypeNames_[gm]
                << ". Available geometry modes include "
                << geometryModeTypeNames_ << exit(FatalError);
        }
    }

    coordSys_ = cylindricalCS("rotorCoordSys", origin, axis, refDir, false);

    const scalar sumArea = gSum(area_);
    const scalar diameter = Foam::sqrt(4.0*sumArea/mathematical::pi);
    Info<< "    Rotor gometry:" << nl
        << "    - disk diameter = " << diameter << nl
        << "    - disk area     = " << sumArea << nl
        << "    - origin        = " << coordSys_.origin() << nl
        << "    - r-axis        = " << coordSys_.R().e1() << nl
        << "    - psi-axis      = " << coordSys_.R().e2() << nl
        << "    - z-axis        = " << coordSys_.R().e3() << endl;
}


void Foam::fv::rotorDiskSource::constructGeometry()
{
    const vectorField& C = mesh_.C();

    forAll(cells_, i)
    {
        if (area_[i] > ROOTVSMALL)
        {
            const label cellI = cells_[i];

            // position in (planar) rotor co-ordinate system
            x_[i] = coordSys_.localPosition(C[cellI]);

            // cache max radius
            rMax_ = max(rMax_, x_[i].x());

            // swept angle relative to rDir axis [radians] in range 0 -> 2*pi
            scalar psi = x_[i].y();

            // blade flap angle [radians]
            scalar beta =
                flap_.beta0 - flap_.beta1c*cos(psi) - flap_.beta2s*sin(psi);

            // determine rotation tensor to convert from planar system into the
            // rotor cone system
            scalar c = cos(beta);
            scalar s = sin(beta);
            R_[i] = tensor(c, 0, -s, 0, 1, 0, s, 0, c);
            invR_[i] = R_[i].T();
        }
    }
}


Foam::tmp<Foam::vectorField> Foam::fv::rotorDiskSource::inflowVelocity
(
    const volVectorField& U
) const
{
    switch (inletFlow_)
    {
        case ifFixed:
        case ifSurfaceNormal:
        {
            return tmp<vectorField>
            (
                new vectorField(mesh_.nCells(), inletVelocity_)
            );

            break;
        }
        case ifLocal:
        {
            return U.internalField();

            break;
        }
        default:
        {
            FatalErrorIn
            (
                "Foam::tmp<Foam::vectorField> "
                "Foam::fv::rotorDiskSource::inflowVelocity"
                "(const volVectorField&) const"
            )   << "Unknown inlet flow specification" << abort(FatalError);
        }
    }

    return tmp<vectorField>(new vectorField(mesh_.nCells(), vector::zero));
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

Foam::fv::rotorDiskSource::rotorDiskSource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh

)
:
    option(name, modelType, dict, mesh),
    rhoRef_(1.0),
    omega_(0.0),
    nBlades_(0),
    inletFlow_(ifLocal),
    inletVelocity_(vector::zero),
    tipEffect_(1.0),
    flap_(),
    x_(cells_.size(), vector::zero),
    R_(cells_.size(), I),
    invR_(cells_.size(), I),
    area_(cells_.size(), 0.0),
    coordSys_(false),
    localAxesRotation_(),
    rMax_(0.0),
    trim_(trimModel::New(*this, coeffs_)),
    blade_(coeffs_.subDict("blade")),
    profiles_(coeffs_.subDict("profiles"))
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fv::rotorDiskSource::~rotorDiskSource()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::rotorDiskSource::addSup
(
    fvMatrix<vector>& eqn,
    const label fieldI
)
{
    volVectorField force
    (
        IOobject
        (
            name_ + ":rotorForce",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedVector
        (
            "zero",
            eqn.dimensions()/dimVolume,
            vector::zero
        )
    );

    // Read the reference density for incompressible flow
    coeffs_.lookup("rhoRef") >> rhoRef_;

    const vectorField Uin(inflowVelocity(eqn.psi()));
    trim_->correct(Uin, force);
    calculate(geometricOneField(), Uin, trim_->thetag(), force);

    // Add source to rhs of eqn
    eqn -= force;

    if (mesh_.time().outputTime())
    {
        force.write();
    }
}


void Foam::fv::rotorDiskSource::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldI
)
{
    volVectorField force
    (
        IOobject
        (
            name_ + ":rotorForce",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedVector
        (
            "zero",
            eqn.dimensions()/dimVolume,
            vector::zero
        )
    );

    const vectorField Uin(inflowVelocity(eqn.psi()));
    trim_->correct(rho, Uin, force);
    calculate(rho, Uin, trim_->thetag(), force);

    // Add source to rhs of eqn
    eqn -= force;

    if (mesh_.time().outputTime())
    {
        force.write();
    }
}


void Foam::fv::rotorDiskSource::writeData(Ostream& os) const
{
    os  << indent << name_ << endl;
    dict_.write(os);
}


bool Foam::fv::rotorDiskSource::read(const dictionary& dict)
{
    if (option::read(dict))
    {
        coeffs_.lookup("fieldNames") >> fieldNames_;
        applied_.setSize(fieldNames_.size(), false);

        // read co-ordinate system/geometry invariant properties
        scalar rpm(readScalar(coeffs_.lookup("rpm")));
        omega_ = rpm/60.0*mathematical::twoPi;

        coeffs_.lookup("nBlades") >> nBlades_;

        inletFlow_ = inletFlowTypeNames_.read(coeffs_.lookup("inletFlowType"));

        coeffs_.lookup("tipEffect") >> tipEffect_;

        const dictionary& flapCoeffs(coeffs_.subDict("flapCoeffs"));
        flapCoeffs.lookup("beta0") >> flap_.beta0;
        flapCoeffs.lookup("beta1c") >> flap_.beta1c;
        flapCoeffs.lookup("beta2s") >> flap_.beta2s;
        flap_.beta0 = degToRad(flap_.beta0);
        flap_.beta1c = degToRad(flap_.beta1c);
        flap_.beta2s = degToRad(flap_.beta2s);


        // create co-ordinate system
        createCoordinateSystem();

        // read co-odinate system dependent properties
        checkData();

        constructGeometry();

        trim_->read(coeffs_);

        if (debug)
        {
            writeField("thetag", trim_->thetag()(), true);
            writeField("faceArea", area_, true);
        }

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
