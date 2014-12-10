/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

Application
    steadyParticleTracks

Description
    Generates a VTK file of particle tracks for cases that were computed using
    a steady-state cloud
    NOTE: case must be re-constructed (if running in parallel) before use

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Cloud.H"
#include "IOdictionary.H"
#include "fvMesh.H"
#include "Time.H"
#include "timeSelector.H"
#include "OFstream.H"
#include "passiveParticleCloud.H"

#include "SortableList.H"
#include "IOobjectList.H"
#include "PtrList.H"
#include "Field.H"
#include "steadyParticleTracksTemplates.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam;

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

label validateFields
(
    const List<word>& userFields,
    const IOobjectList& cloudObjs
)
{
    List<bool> ok(userFields.size(), false);

    forAll(userFields, i)
    {
        ok[i] = ok[i] || fieldOk<label>(cloudObjs, userFields[i]);
        ok[i] = ok[i] || fieldOk<scalar>(cloudObjs, userFields[i]);
        ok[i] = ok[i] || fieldOk<vector>(cloudObjs, userFields[i]);
        ok[i] = ok[i] || fieldOk<sphericalTensor>(cloudObjs, userFields[i]);
        ok[i] = ok[i] || fieldOk<symmTensor>(cloudObjs, userFields[i]);
        ok[i] = ok[i] || fieldOk<tensor>(cloudObjs, userFields[i]);
    }

    label nOk = 0;
    forAll(ok, i)
    {
        if (ok[i])
        {
            nOk++;
        }
        else
        {
            Info << "\n*** Warning: user specified field '" << userFields[i]
                 << "' unavailable" << endl;
        }
    }

    return nOk;
}


template<>
void writeVTK(OFstream& os, const label& value)
{
    os  << value;
}


template<>
void writeVTK(OFstream& os, const scalar& value)
{
    os  << value;
}

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noParallel();
    timeSelector::addOptions();
    #include "addRegionOption.H"
    #include "addDictOption.H"

    #include "setRootCase.H"

    #include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
    #include "createNamedMesh.H"
    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    fileName vtkPath(runTime.path()/"VTK");
    mkDir(vtkPath);

    typedef HashTable<label, labelPair, labelPair::Hash<> > trackTableType;

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.timeName() << endl;

        fileName vtkTimePath(runTime.path()/"VTK"/runTime.timeName());
        mkDir(vtkTimePath);

        Info<< "    Reading particle positions" << endl;

        PtrList<passiveParticle> particles(0);

        // transfer particles to (more convenient) list
        {
            passiveParticleCloud ppc(mesh, cloudName);
            Info<< "\n    Read " << returnReduce(ppc.size(), sumOp<label>())
                << " particles" << endl;

            particles.setSize(ppc.size());

            label i = 0;
            forAllIter(passiveParticleCloud, ppc, iter)
            {
                particles.set(i++, ppc.remove(&iter()));
            }

            // myCloud should now be empty
        }

        List<label> particleToTrack(particles.size());
        label nTracks = 0;

        {
            trackTableType trackTable;
            forAll(particles, i)
            {
                const label origProc = particles[i].origProc();
                const label origId = particles[i].origId();

                const trackTableType::const_iterator& iter =
                    trackTable.find(labelPair(origProc, origId));

                if (iter == trackTable.end())
                {
                    particleToTrack[i] = nTracks;
                    trackTable.insert(labelPair(origProc, origId), nTracks);
                    nTracks++;
                }
                else
                {
                    particleToTrack[i] = iter();
                }
            }
        }


        if (nTracks == 0)
        {
            Info<< "\n    No track data" << endl;
        }
        else
        {
            Info<< "\n    Generating " << nTracks << " tracks" << endl;

            // determine length of each track
            labelList trackLengths(nTracks, 0);
            forAll(particleToTrack, i)
            {
                const label trackI = particleToTrack[i];
                trackLengths[trackI]++;
            }

            // particle "age" property used to sort the tracks
            List<SortableList<scalar> > agePerTrack(nTracks);
            List<List<label> > particleMap(nTracks);

            forAll(trackLengths, i)
            {
                const label length = trackLengths[i];
                agePerTrack[i].setSize(length);
                particleMap[i].setSize(length);
            }

            // store the particle age per track
            IOobjectList cloudObjs
            (
                mesh,
                runTime.timeName(),
                cloud::prefix/cloudName
            );

            // TODO: gather age across all procs
            {
                tmp<scalarField> tage =
                    readParticleField<scalar>("age", cloudObjs);
                const scalarField& age = tage();
                List<label> trackSamples(nTracks, 0);
                forAll(particleToTrack, i)
                {
                    const label trackI = particleToTrack[i];
                    const label sampleI = trackSamples[trackI];
                    agePerTrack[trackI][sampleI] = age[i];
                    particleMap[trackI][sampleI] = i;
                    trackSamples[trackI]++;
                }
                tage.clear();
            }


            if (Pstream::master())
            {
                OFstream os(vtkTimePath/"particleTracks.vtk");

                Info<< "\n    Writing particle tracks to " << os.name() << endl;

                label nPoints = sum(trackLengths);

                os  << "# vtk DataFile Version 2.0" << nl
                    << "particleTracks" << nl
                    << "ASCII" << nl
                    << "DATASET POLYDATA" << nl
                    << "POINTS " << nPoints << " float" << nl;

                Info<< "\n    Writing points" << endl;

                {
                    forAll(agePerTrack, i)
                    {
                        agePerTrack[i].sort();

                        const labelList& ids = agePerTrack[i].indices();
                        labelList& particleIds = particleMap[i];

                        {
                            // update addressing
                            List<label> sortedIds(ids);
                            forAll(sortedIds, j)
                            {
                                sortedIds[j] = particleIds[ids[j]];
                            }
                            particleIds = sortedIds;
                        }

                        forAll(ids, j)
                        {
                            const label localId = particleIds[j];
                            const vector& pos = particles[localId].position();
                            os  << pos.x() << ' ' << pos.y() << ' ' << pos.z()
                                << nl;
                        }
                    }
                }


                // write track (line) connectivity to file

                Info<< "\n    Writing track lines" << endl;
                os  << "\nLINES " << nTracks << ' ' << nPoints + nTracks << nl;

                // Write ids of track points to file
                {
                    label globalPtI = 0;
                    forAll(particleMap, i)
                    {
                        os  << particleMap[i].size() << nl;

                        forAll(particleMap[i], j)
                        {
                            os  << ' ' << globalPtI++;

                            if (((j + 1) % 10 == 0) && (j != 0))
                            {
                                os  << nl;
                            }
                        }

                        os  << nl;
                    }
                }


                const label nFields = validateFields(userFields, cloudObjs);

                os  << "POINT_DATA " << nPoints << nl
                    << "FIELD attributes " << nFields << nl;

                Info<< "\n    Processing fields" << nl << endl;

                processFields<label>(os, particleMap, userFields, cloudObjs);
                processFields<scalar>(os, particleMap, userFields, cloudObjs);
                processFields<vector>(os, particleMap, userFields, cloudObjs);
                processFields<sphericalTensor>
                    (os, particleMap, userFields, cloudObjs);
                processFields<symmTensor>
                    (os, particleMap, userFields, cloudObjs);
                processFields<tensor>(os, particleMap, userFields, cloudObjs);

            }
        }
        Info<< endl;
    }

    Info<< "\ndone" << endl;

    return 0;
}


// ************************************************************************* //
