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

#include "rotorDiskSource.H"
#include "trimModel.H"
#include "fvMatrices.H"
#include "geometricOneField.H"
#include "syncTools.H"
#include "axesRotation.H"
#include "addToRunTimeSelectionTable.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
    namespace fv
    {
        defineTypeNameAndDebug(rotorDiskSource, 0);
        addToRunTimeSelectionTable(fvModel, rotorDiskSource, dictionary);
    }

    template<> const char*
    NamedEnum<fv::rotorDiskSource::geometryModeType, 2>::names[] =
    {
        "auto",
        "specified"
    };

    const NamedEnum<fv::rotorDiskSource::geometryModeType, 2>
        fv::rotorDiskSource::geometryModeTypeNames_;

    template<> const char*
    NamedEnum<fv::rotorDiskSource::inletFlowType, 3>::names[] =
    {
        "fixed",
        "surfaceNormal",
        "local"
    };

    const NamedEnum<fv::rotorDiskSource::inletFlowType, 3>
        fv::rotorDiskSource::inletFlowTypeNames_;
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::fv::rotorDiskSource::readCoeffs()
{
    UName_ = coeffs().lookupOrDefault<word>("U", "U");

    // Read co-ordinate system/geometry invariant properties
    scalar rpm(coeffs().lookup<scalar>("rpm"));
    omega_ = rpm/60.0*mathematical::twoPi;

    coeffs().lookup("nBlades") >> nBlades_;

    inletFlow_ = inletFlowTypeNames_.read(coeffs().lookup("inletFlowType"));

    coeffs().lookup("tipEffect") >> tipEffect_;

    const dictionary& flapCoeffs(coeffs().subDict("flapCoeffs"));
    flapCoeffs.lookup("beta0") >> flap_.beta0;
    flapCoeffs.lookup("beta1c") >> flap_.beta1c;
    flapCoeffs.lookup("beta2s") >> flap_.beta2s;
    flap_.beta0 = degToRad(flap_.beta0);
    flap_.beta1c = degToRad(flap_.beta1c);
    flap_.beta2s = degToRad(flap_.beta2s);

    // Create co-ordinate system
    createCoordinateSystem();

    // Read co-odinate system dependent properties
    checkData();

    constructGeometry();

    trim_->read(coeffs());

    if (debug)
    {
        writeField("thetag", trim_->thetag()(), true);
        writeField("faceArea", area_, true);
    }
}


void Foam::fv::rotorDiskSource::checkData()
{
    // Set inflow type
    switch (set_.selectionMode())
    {
        case fvCellSet::selectionModeType::cellSet:
        case fvCellSet::selectionModeType::cellZone:
        case fvCellSet::selectionModeType::all:
        {
            // Set the profile ID for each blade section
            profiles_.connectBlades(blade_.profileName(), blade_.profileID());
            switch (inletFlow_)
            {
                case inletFlowType::fixed:
                {
                    coeffs().lookup("inletVelocity") >> inletVelocity_;
                    break;
                }
                case inletFlowType::surfaceNormal:
                {
                    scalar UIn
                    (
                        coeffs().lookup<scalar>("inletNormalVelocity")
                    );
                    inletVelocity_ = -coordSys_.R().e3()*UIn;
                    break;
                }
                case inletFlowType::local:
                {
                    break;
                }
                default:
                {
                    FatalErrorInFunction
                        << "Unknown inlet velocity type" << abort(FatalError);
                }
            }
            break;
        }
        default:
        {
            FatalErrorInFunction
                << "Source cannot be used with '"
                << fvCellSet::selectionModeTypeNames_[set_.selectionMode()]
                << "' mode.  Please use one of: " << nl
                << fvCellSet::selectionModeTypeNames_
                   [fvCellSet::selectionModeType::cellSet] << nl
                << fvCellSet::selectionModeTypeNames_
                   [fvCellSet::selectionModeType::cellZone] << nl
                << fvCellSet::selectionModeTypeNames_
                   [fvCellSet::selectionModeType::all]
                << exit(FatalError);
        }
    }
}


void Foam::fv::rotorDiskSource::setFaceArea(vector& axis, const bool correct)
{
    area_ = 0.0;

    static const scalar tol = 0.8;

    const label nInternalFaces = mesh().nInternalFaces();
    const polyBoundaryMesh& pbm = mesh().boundaryMesh();
    const vectorField& Sf = mesh().Sf();
    const scalarField& magSf = mesh().magSf();

    vector n = Zero;

    // Calculate cell addressing for selected cells
    labelList cellAddr(mesh().nCells(), -1);
    UIndirectList<label>(cellAddr, set_.cells()) =
        identity(set_.cells().size());
    labelList nbrFaceCellAddr(mesh().nFaces() - nInternalFaces, -1);
    forAll(pbm, patchi)
    {
        const polyPatch& pp = pbm[patchi];

        if (pp.coupled())
        {
            forAll(pp, i)
            {
                label facei = pp.start() + i;
                label nbrFacei = facei - nInternalFaces;
                label own = mesh().faceOwner()[facei];
                nbrFaceCellAddr[nbrFacei] = cellAddr[own];
            }
        }
    }

    // Correct for parallel running
    syncTools::swapBoundaryFaceList(mesh(), nbrFaceCellAddr);

    // Add internal field contributions
    for (label facei = 0; facei < nInternalFaces; facei++)
    {
        const label own = cellAddr[mesh().faceOwner()[facei]];
        const label nbr = cellAddr[mesh().faceNeighbour()[facei]];

        if ((own != -1) && (nbr == -1))
        {
            vector nf = Sf[facei]/magSf[facei];

            if ((nf & axis) > tol)
            {
                area_[own] += magSf[facei];
                n += Sf[facei];
            }
        }
        else if ((own == -1) && (nbr != -1))
        {
            vector nf = Sf[facei]/magSf[facei];

            if ((-nf & axis) > tol)
            {
                area_[nbr] += magSf[facei];
                n -= Sf[facei];
            }
        }
    }


    // Add boundary contributions
    forAll(pbm, patchi)
    {
        const polyPatch& pp = pbm[patchi];
        const vectorField& Sfp = mesh().Sf().boundaryField()[patchi];
        const scalarField& magSfp = mesh().magSf().boundaryField()[patchi];

        if (pp.coupled())
        {
            forAll(pp, j)
            {
                const label facei = pp.start() + j;
                const label own = cellAddr[mesh().faceOwner()[facei]];
                const label nbr = nbrFaceCellAddr[facei - nInternalFaces];
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
                const label facei = pp.start() + j;
                const label own = cellAddr[mesh().faceOwner()[facei]];
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
                name() + ":area",
                mesh().time().name(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedScalar(dimArea, 0)
        );
        UIndirectList<scalar>(area.primitiveField(), set_.cells()) = area_;

        Info<< type() << ": " << name() << " writing field " << area.name()
            << endl;

        area.write();
    }
}


void Foam::fv::rotorDiskSource::createCoordinateSystem()
{
    // Construct the local rotor co-prdinate system
    vector origin(Zero);
    vector axis(Zero);
    vector refDir(Zero);

    geometryModeType gm =
        geometryModeTypeNames_.read(coeffs().lookup("geometryMode"));

    switch (gm)
    {
        case geometryModeType::automatic:
        {
            // Determine rotation origin (cell volume weighted)
            scalar sumV = 0.0;
            const scalarField& V = mesh().V();
            const vectorField& C = mesh().C();

            const labelList& cells = set_.cells();

            forAll(cells, i)
            {
                const label celli = cells[i];
                sumV += V[celli];
                origin += V[celli]*C[celli];
            }
            reduce(origin, sumOp<vector>());
            reduce(sumV, sumOp<scalar>());
            origin /= sumV;

            // Determine first radial vector
            vector dx1(Zero);
            scalar magR = -great;
            forAll(cells, i)
            {
                const label celli = cells[i];
                vector test = C[celli] - origin;
                if (mag(test) > magR)
                {
                    dx1 = test;
                    magR = mag(test);
                }
            }
            reduce(dx1, maxMagSqrOp<vector>());
            magR = mag(dx1);

            // Determine second radial vector and cross to determine axis
            forAll(cells, i)
            {
                const label celli = cells[i];
                vector dx2 = C[celli] - origin;
                if (mag(dx2) > 0.5*magR)
                {
                    axis = dx1 ^ dx2;
                    if (mag(axis) > small)
                    {
                        break;
                    }
                }
            }
            reduce(axis, maxMagSqrOp<vector>());
            axis /= mag(axis);

            // Correct the axis direction using a point above the rotor
            {
                vector pointAbove(coeffs().lookup("pointAbove"));
                vector dir = pointAbove - origin;
                dir /= mag(dir);
                if ((dir & axis) < 0)
                {
                    axis *= -1.0;
                }
            }

            coeffs().lookup("refDirection") >> refDir;

            cylindrical_.reset
            (
                new cylindrical(axis, origin, UIndirectList<vector>(C, cells)())
            );

            // Set the face areas and apply correction to calculated axis
            // e.g. if cellZone is more than a single layer in thickness
            setFaceArea(axis, true);

            break;
        }
        case geometryModeType::specified:
        {
            coeffs().lookup("origin") >> origin;
            coeffs().lookup("axis") >> axis;
            coeffs().lookup("refDirection") >> refDir;

            cylindrical_.reset
            (
                new cylindrical
                (
                    axis,
                    origin,
                    UIndirectList<vector>(mesh().C(), set_.cells())()
                )
            );

            setFaceArea(axis, false);

            break;
        }
    }

    coordSys_ = coordinateSystems::cylindrical
    (
        "rotorCoordSys",
        origin,
        axis,
        refDir,
        false
    );

    const scalar sumArea = gSum(area_);
    const scalar diameter = Foam::sqrt(4.0*sumArea/mathematical::pi);
    Info<< "    Rotor geometry:" << nl
        << "    - disk diameter = " << diameter << nl
        << "    - disk area     = " << sumArea << nl
        << "    - origin        = " << coordSys_.origin() << nl
        << "    - r-axis        = " << coordSys_.R().e1() << nl
        << "    - psi-axis      = " << coordSys_.R().e2() << nl
        << "    - z-axis        = " << coordSys_.R().e3() << endl;
}


void Foam::fv::rotorDiskSource::constructGeometry()
{
    const vectorField& C = mesh().C();

    const labelList& cells = set_.cells();

    forAll(cells, i)
    {
        if (area_[i] > rootVSmall)
        {
            const label celli = cells[i];

            // Position in (planar) rotor co-ordinate system
            x_[i] = coordSys_.localPosition(C[celli]);

            // Cache max radius
            rMax_ = max(rMax_, x_[i].x());

            // Swept angle relative to rDir axis [radians] in range 0 -> 2*pi
            scalar psi = x_[i].y();

            // Blade flap angle [radians]
            scalar beta =
                flap_.beta0 - flap_.beta1c*cos(psi) - flap_.beta2s*sin(psi);

            // Determine rotation tensor to convert from planar system into the
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
        case inletFlowType::fixed:
        case inletFlowType::surfaceNormal:
        {
            return tmp<vectorField>
            (
                new vectorField(mesh().nCells(), inletVelocity_)
            );

            break;
        }
        case inletFlowType::local:
        {
            return U.primitiveField();

            break;
        }
        default:
        {
            FatalErrorInFunction
                << "Unknown inlet flow specification" << abort(FatalError);
        }
    }

    return tmp<vectorField>(new vectorField(mesh().nCells(), Zero));
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
    fvModel(name, modelType, dict, mesh),
    set_(mesh, coeffs()),
    UName_(word::null),
    omega_(0),
    nBlades_(0),
    inletFlow_(inletFlowType::local),
    inletVelocity_(Zero),
    tipEffect_(1),
    flap_(),
    x_(set_.cells().size(), Zero),
    R_(set_.cells().size(), I),
    invR_(set_.cells().size(), I),
    area_(set_.cells().size(), Zero),
    coordSys_("rotorCoordSys", vector::zero, axesRotation(sphericalTensor::I)),
    cylindrical_(),
    rMax_(0),
    trim_(trimModel::New(*this, coeffs())),
    blade_(coeffs().subDict("blade")),
    profiles_(coeffs().subDict("profiles")),
    rhoRef_(1)
{
    readCoeffs();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fv::rotorDiskSource::~rotorDiskSource()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::wordList Foam::fv::rotorDiskSource::addSupFields() const
{
    return wordList(1, UName_);
}


void Foam::fv::rotorDiskSource::addSup
(
    fvMatrix<vector>& eqn,
    const word& fieldName
) const
{
    volVectorField force
    (
        IOobject
        (
            name() + ":rotorForce",
            mesh().time().name(),
            mesh()
        ),
        mesh(),
        dimensionedVector
        (
            "zero",
            eqn.dimensions()/dimVolume,
            Zero
        )
    );

    // Read the reference density for incompressible flow
    coeffs().lookup("rhoRef") >> rhoRef_;

    const vectorField Uin(inflowVelocity(eqn.psi()));
    trim_->correct(Uin, force);
    calculate(geometricOneField(), Uin, trim_->thetag(), force);

    // Add source to rhs of eqn
    eqn -= force;

    if (mesh().time().writeTime())
    {
        force.write();
    }
}


void Foam::fv::rotorDiskSource::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const word& fieldName
) const
{
    volVectorField force
    (
        IOobject
        (
            name() + ":rotorForce",
            mesh().time().name(),
            mesh()
        ),
        mesh(),
        dimensionedVector
        (
            "zero",
            eqn.dimensions()/dimVolume,
            Zero
        )
    );

    const vectorField Uin(inflowVelocity(eqn.psi()));
    trim_->correct(rho, Uin, force);
    calculate(rho, Uin, trim_->thetag(), force);

    // Add source to rhs of eqn
    eqn -= force;

    if (mesh().time().writeTime())
    {
        force.write();
    }
}


bool Foam::fv::rotorDiskSource::movePoints()
{
    set_.movePoints();
    return true;
}


void Foam::fv::rotorDiskSource::topoChange(const polyTopoChangeMap& map)
{
    set_.topoChange(map);
}


void Foam::fv::rotorDiskSource::mapMesh(const polyMeshMap& map)
{
    set_.mapMesh(map);
}


void Foam::fv::rotorDiskSource::distribute(const polyDistributionMap& map)
{
    set_.distribute(map);
}


bool Foam::fv::rotorDiskSource::read(const dictionary& dict)
{
    if (fvModel::read(dict))
    {
        set_.read(coeffs());
        readCoeffs();
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
