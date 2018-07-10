/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "MapLagrangianFields.H"
#include "passiveParticleCloud.H"
#include "meshSearch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

static const scalar perturbFactor = 1e-6;


// Special version of findCell that generates a cell guaranteed to be
// compatible with tracking.
static label findCell(const Cloud<passiveParticle>& cloud, const point& pt)
{
    label celli = -1;
    label tetFacei = -1;
    label tetPtI = -1;

    const polyMesh& mesh = cloud.pMesh();

    mesh.findCellFacePt(pt, celli, tetFacei, tetPtI);

    if (celli >= 0)
    {
        return celli;
    }
    else
    {
        // See if particle on face by finding nearest face and shifting
        // particle.

        meshSearch meshSearcher
        (
            mesh,
            polyMesh::FACE_PLANES    // no decomposition needed
        );

        label facei = meshSearcher.findNearestBoundaryFace(pt);

        if (facei >= 0)
        {
            const point& cc = mesh.cellCentres()[mesh.faceOwner()[facei]];

            const point perturbPt = (1-perturbFactor)*pt+perturbFactor*cc;

            mesh.findCellFacePt(perturbPt, celli, tetFacei, tetPtI);

            return celli;
        }
    }

    return -1;
}


void mapLagrangian(const meshToMesh0& meshToMesh0Interp)
{
    // Determine which particles are in meshTarget
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // target to source cell map
    const labelList& cellAddressing = meshToMesh0Interp.cellAddressing();

    // Invert celladdressing to get source to target(s).
    // Note: could use sparse addressing but that is too storage inefficient
    // (Map<labelList>)
    labelListList sourceToTargets
    (
        invertOneToMany(meshToMesh0Interp.fromMesh().nCells(), cellAddressing)
    );

    const fvMesh& meshSource = meshToMesh0Interp.fromMesh();
    const fvMesh& meshTarget = meshToMesh0Interp.toMesh();


    fileNameList cloudDirs
    (
        readDir
        (
            meshSource.time().timePath()/cloud::prefix,
            fileName::DIRECTORY
        )
    );

    forAll(cloudDirs, cloudI)
    {
        // Search for list of lagrangian objects for this time
        IOobjectList objects
        (
            meshSource,
            meshSource.time().timeName(),
            cloud::prefix/cloudDirs[cloudI]
        );

        IOobject* positionsPtr = objects.lookup("positions");

        if (positionsPtr)
        {
            Info<< nl << "    processing cloud " << cloudDirs[cloudI] << endl;

            // Read positions & cell
            passiveParticleCloud sourceParcels
            (
                meshSource,
                cloudDirs[cloudI],
                false
            );
            Info<< "    read " << sourceParcels.size()
                << " parcels from source mesh." << endl;

            // Construct empty target cloud
            passiveParticleCloud targetParcels
            (
                meshTarget,
                cloudDirs[cloudI],
                IDLList<passiveParticle>()
            );

            passiveParticle::trackingData td(targetParcels);

            label sourceParticleI = 0;

            // Indices of source particles that get added to targetParcels
            DynamicList<label> addParticles(sourceParcels.size());

            // Unmapped particles
            labelHashSet unmappedSource(sourceParcels.size());


            // Initial: track from fine-mesh cell centre to particle position
            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            // This requires there to be no boundary in the way.


            forAllConstIter(Cloud<passiveParticle>, sourceParcels, iter)
            {
                bool foundCell = false;

                // Assume that cell from read parcel is the correct one...
                if (iter().cell() >= 0)
                {
                    const labelList& targetCells =
                        sourceToTargets[iter().cell()];

                    // Particle probably in one of the targetcells. Try
                    // all by tracking from their cell centre to the parcel
                    // position.

                    forAll(targetCells, i)
                    {
                        // Track from its cellcentre to position to make sure.
                        autoPtr<passiveParticle> newPtr
                        (
                            new passiveParticle
                            (
                                meshTarget,
                                barycentric(1, 0, 0, 0),
                                targetCells[i],
                                meshTarget.cells()[targetCells[i]][0],
                                1
                            )
                        );
                        passiveParticle& newP = newPtr();

                        newP.track(iter().position() - newP.position(), 0);

                        if (!newP.onFace())
                        {
                            // Hit position.
                            foundCell = true;
                            addParticles.append(sourceParticleI);
                            targetParcels.addParticle(newPtr.ptr());
                            break;
                        }
                    }
                }

                if (!foundCell)
                {
                    // Store for closer analysis
                    unmappedSource.insert(sourceParticleI);
                }

                sourceParticleI++;
            }

            Info<< "    after meshToMesh0 addressing found "
                << targetParcels.size()
                << " parcels in target mesh." << endl;


            // Do closer inspection for unmapped particles
            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            if (unmappedSource.size())
            {
                sourceParticleI = 0;

                forAllIter(Cloud<passiveParticle>, sourceParcels, iter)
                {
                    if (unmappedSource.found(sourceParticleI))
                    {
                        label targetCell =
                            findCell(targetParcels, iter().position());

                        if (targetCell >= 0)
                        {
                            unmappedSource.erase(sourceParticleI);
                            addParticles.append(sourceParticleI);
                            targetParcels.addParticle
                            (
                                new passiveParticle
                                (
                                    meshTarget,
                                    iter().position(),
                                    targetCell
                                )
                            );
                            sourceParcels.remove(&iter());
                        }
                    }
                    sourceParticleI++;
                }
            }
            addParticles.shrink();

            Info<< "    after additional mesh searching found "
                << targetParcels.size() << " parcels in target mesh." << endl;

            if (addParticles.size())
            {
                IOPosition<passiveParticleCloud>(targetParcels).write();

                // addParticles now contains the indices of the sourceMesh
                // particles that were appended to the target mesh.

                // Map lagrangian fields
                // ~~~~~~~~~~~~~~~~~~~~~

                MapLagrangianFields<label>
                (cloudDirs[cloudI], objects, meshToMesh0Interp, addParticles);
                MapLagrangianFields<scalar>
                (cloudDirs[cloudI], objects, meshToMesh0Interp, addParticles);
                MapLagrangianFields<vector>
                (cloudDirs[cloudI], objects, meshToMesh0Interp, addParticles);
                MapLagrangianFields<sphericalTensor>
                (cloudDirs[cloudI], objects, meshToMesh0Interp, addParticles);
                MapLagrangianFields<symmTensor>
                (cloudDirs[cloudI], objects, meshToMesh0Interp, addParticles);
                MapLagrangianFields<tensor>
                (cloudDirs[cloudI], objects, meshToMesh0Interp, addParticles);
            }
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
