/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2024 OpenFOAM Foundation
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

#include "addToRunTimeSelectionTable.H"
#include "powerLawLopesdaCosta.H"
#include "geometricOneField.H"
#include "fvMatrices.H"
#include "triSurfaceMesh.H"
#include "triSurfaceTools.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace porosityModels
    {
        defineTypeNameAndDebug(powerLawLopesdaCosta, 0);
        addToRunTimeSelectionTable(porosityModel, powerLawLopesdaCosta, mesh);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::porosityModels::powerLawLopesdaCostaZone::powerLawLopesdaCostaZone
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict,
    const dictionary& coeffDict
)
:
    zoneName_(name + ":porousZone")
{
    // Vertical direction within porous region
    vector zDir(coeffDict.lookup("zDir"));

    // Span of the search for the cells in the porous region
    scalar searchSpan(coeffDict.lookup<scalar>("searchSpan"));

    // Top surface file name defining the extent of the porous region
    word topSurfaceFileName(coeffDict.lookup("topSurface"));

    // List of the ground patches defining the lower surface
    // of the porous region
    List<wordRe> groundPatches(coeffDict.lookup("groundPatches"));

    // Functional form of the porosity surface area per unit volume as a
    // function of the normalised vertical position
    autoPtr<Function1<scalar>> AvFunc
    (
        Function1<scalar>::New
        (
            dict.lookupEntryBackwardsCompatible
            (
                {"Av", "Sigma"},
                false,
                true
            ).keyword(),
            dimLength,
            dimArea/dimVolume,
            dict
        )
    );

    // Searchable triSurface for the top of the porous region
    triSurfaceMesh searchSurf
    (
        IOobject
        (
            topSurfaceFileName,
            mesh.time().constant(),
            searchableSurface::geometryDir(mesh.time()),
            mesh.time()
        )
    );

    // Limit the porous cell search to those within the lateral and vertical
    // extent of the top surface

    boundBox surfBounds(searchSurf.points());
    boundBox searchBounds
    (
        surfBounds.min() - searchSpan*zDir, surfBounds.max()
    );

    const pointField& C = mesh.cellCentres();

    // List of cells within the porous region
    labelList porousCells(C.size());
    label porousCelli = 0;

    forAll(C, celli)
    {
        if (searchBounds.contains(C[celli]))
        {
            porousCells[porousCelli++] = celli;
        }
    }

    porousCells.setSize(porousCelli);

    pointField start(porousCelli);
    pointField end(porousCelli);

    forAll(porousCells, porousCelli)
    {
        start[porousCelli] = C[porousCells[porousCelli]];
        end[porousCelli] = start[porousCelli] + searchSpan*zDir;
    }

    // Field of distance between the cell centre and the corresponding point
    // on the porous region top surface
    scalarField zTop(porousCelli);

    // Search the porous cells for those with a corresponding point on the
    // porous region top surface
    List<pointIndexHit> hitInfo;
    searchSurf.findLine(start, end, hitInfo);

    porousCelli = 0;

    forAll(porousCells, celli)
    {
        const pointIndexHit& hit = hitInfo[celli];

        if (hit.hit())
        {
            porousCells[porousCelli] = porousCells[celli];

            zTop[porousCelli] =
                (hit.hitPoint() - C[porousCells[porousCelli]]) & zDir;

            porousCelli++;
        }
    }

    // Resize the porous cells list and height field
    porousCells.setSize(porousCelli);
    zTop.setSize(porousCelli);

    // Create a triSurface for the ground patch(s)
    triSurface groundSurface
    (
        triSurfaceTools::triangulate
        (
            mesh.boundaryMesh(),
            mesh.boundaryMesh().patchSet(groundPatches),
            searchBounds
        )
    );

    // Combine the ground triSurfaces across all processors
    if (Pstream::parRun())
    {
        List<List<labelledTri>> groundSurfaceProcTris(Pstream::nProcs());
        List<pointField> groundSurfaceProcPoints(Pstream::nProcs());

        groundSurfaceProcTris[Pstream::myProcNo()] = groundSurface;
        groundSurfaceProcPoints[Pstream::myProcNo()] = groundSurface.points();

        Pstream::gatherList(groundSurfaceProcTris);
        Pstream::scatterList(groundSurfaceProcTris);
        Pstream::gatherList(groundSurfaceProcPoints);
        Pstream::scatterList(groundSurfaceProcPoints);

        label nTris = 0;
        forAll(groundSurfaceProcTris, i)
        {
            nTris += groundSurfaceProcTris[i].size();
        }

        List<labelledTri> groundSurfaceTris(nTris);
        label trii = 0;
        label offset = 0;
        forAll(groundSurfaceProcTris, i)
        {
            forAll(groundSurfaceProcTris[i], j)
            {
                groundSurfaceTris[trii] = groundSurfaceProcTris[i][j];
                groundSurfaceTris[trii][0] += offset;
                groundSurfaceTris[trii][1] += offset;
                groundSurfaceTris[trii][2] += offset;
                trii++;
            }
            offset += groundSurfaceProcPoints[i].size();
        }

        label nPoints = 0;
        forAll(groundSurfaceProcPoints, i)
        {
            nPoints += groundSurfaceProcPoints[i].size();
        }

        pointField groundSurfacePoints(nPoints);

        label pointi = 0;
        forAll(groundSurfaceProcPoints, i)
        {
            forAll(groundSurfaceProcPoints[i], j)
            {
                groundSurfacePoints[pointi++] = groundSurfaceProcPoints[i][j];
            }
        }

        groundSurface = triSurface(groundSurfaceTris, groundSurfacePoints);
    }

    // Create a searchable triSurface for the ground
    triSurfaceSearch groundSearch(groundSurface);

    // Search the porous cells for the corresponding point on the ground

    start.setSize(porousCelli);
    end.setSize(porousCelli);

    forAll(porousCells, porousCelli)
    {
        start[porousCelli] = C[porousCells[porousCelli]];
        end[porousCelli] = start[porousCelli] - searchSpan*zDir;
    }

    groundSearch.findLine(start, end, hitInfo);

    scalarField zBottom(porousCelli);

    forAll(porousCells, porousCelli)
    {
        const pointIndexHit& hit = hitInfo[porousCelli];

        if (hit.hit())
        {
            zBottom[porousCelli] =
                (C[porousCells[porousCelli]] - hit.hitPoint()) & zDir;
        }
    }

    // Create the normalised height field
    const scalarField zNorm(zBottom/(zBottom + zTop));

    // Create the porosity surface area per unit volume zone field
    Av_ = AvFunc->value(zNorm);

    // Create the porous region cellZone and add to the mesh cellZones
    if (!mesh.cellZones().found(zoneName_))
    {
        mesh.cellZones().append
        (
            new cellZone
            (
                zoneName_,
                porousCells,
                mesh.cellZones()
            )
        );
    }
    else
    {
        FatalErrorInFunction
            << "Unable to create porous cellZone " << zoneName_
            << ": zone already exists"
            << abort(FatalError);
    }
}


Foam::porosityModels::powerLawLopesdaCosta::powerLawLopesdaCosta
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict,
    const dictionary& coeffDict,
    const word& dummyCellZoneName
)
:
    powerLawLopesdaCostaZone(name, mesh, dict, coeffDict),
    porosityModel
    (
        name,
        mesh,
        dict,
        coeffDict,
        powerLawLopesdaCostaZone::zoneName_
    ),
    coeffDict_(coeffDict),
    Cd_(coeffDict.lookup<scalar>("Cd")),
    C1_(coeffDict.lookup<scalar>("C1")),
    rhoName_(coeffDict.lookupOrDefault<word>("rho", "rho"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::porosityModels::powerLawLopesdaCosta::~powerLawLopesdaCosta()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::scalarField&
Foam::porosityModels::powerLawLopesdaCostaZone::Av() const
{
    return Av_;
}


void Foam::porosityModels::powerLawLopesdaCosta::calcTransformModelData()
{}


void Foam::porosityModels::powerLawLopesdaCosta::calcForce
(
    const volVectorField& U,
    const volScalarField& rho,
    const volScalarField& mu,
    vectorField& force
) const
{
    scalarField Udiag(U.size(), 0.0);
    const scalarField& V = mesh_.V();

    apply(Udiag, V, rho, U);

    force = Udiag*U;
}


void Foam::porosityModels::powerLawLopesdaCosta::correct
(
    fvVectorMatrix& UEqn
) const
{
    const vectorField& U = UEqn.psi();
    const scalarField& V = mesh_.V();
    scalarField& Udiag = UEqn.diag();

    if (UEqn.dimensions() == dimForce)
    {
        const volScalarField& rho =
            mesh_.lookupObject<volScalarField>(rhoName_);

        apply(Udiag, V, rho, U);
    }
    else
    {
        apply(Udiag, V, geometricOneField(), U);
    }
}


void Foam::porosityModels::powerLawLopesdaCosta::correct
(
    fvVectorMatrix& UEqn,
    const volScalarField& rho,
    const volScalarField& mu
) const
{
    const vectorField& U = UEqn.psi();
    const scalarField& V = mesh_.V();
    scalarField& Udiag = UEqn.diag();

    apply(Udiag, V, rho, U);
}


void Foam::porosityModels::powerLawLopesdaCosta::correct
(
    const fvVectorMatrix& UEqn,
    volTensorField& AU
) const
{
    const vectorField& U = UEqn.psi();

    if (UEqn.dimensions() == dimForce)
    {
        const volScalarField& rho =
            mesh_.lookupObject<volScalarField>(rhoName_);

        apply(AU, rho, U);
    }
    else
    {
        apply(AU, geometricOneField(), U);
    }
}


// ************************************************************************* //
