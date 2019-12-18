/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2019 OpenFOAM Foundation
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

#include "coupleGroupIdentifier.H"
#include "polyMesh.H"
#include "Time.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::label Foam::coupleGroupIdentifier::findOtherPatchID
(
    const polyMesh& mesh,
    const polyPatch& thisPatch
) const
{
    const polyBoundaryMesh& pbm = mesh.boundaryMesh();

    if (!valid())
    {
        FatalErrorInFunction
            << "Invalid coupleGroup patch group"
            << " on patch " << thisPatch.name()
            << " in region " << pbm.mesh().name()
            << exit(FatalError);
    }

    HashTable<labelList, word>::const_iterator fnd =
        pbm.groupPatchIDs().find(name());

    if (fnd == pbm.groupPatchIDs().end())
    {
        if (&mesh == &thisPatch.boundaryMesh().mesh())
        {
            // thisPatch should be in patchGroup
            FatalErrorInFunction
                << "Patch " << thisPatch.name()
                << " should be in patchGroup " << name()
                << " in region " << pbm.mesh().name()
                << exit(FatalError);
        }

        return -1;
    }

    // Mesh has patch group
    const labelList& patchIDs = fnd();

    if (&mesh == &thisPatch.boundaryMesh().mesh())
    {
        if (patchIDs.size() > 2 || patchIDs.size() == 0)
        {
            FatalErrorInFunction
                << "Couple patchGroup " << name()
                << " with contents " << patchIDs
                << " not of size < 2"
                << " on patch " << thisPatch.name()
                << " region " << thisPatch.boundaryMesh().mesh().name()
                << exit(FatalError);

            return -1;
        }

        label index = findIndex(patchIDs, thisPatch.index());

        if (index == -1)
        {
            FatalErrorInFunction
                << "Couple patchGroup " << name()
                << " with contents " << patchIDs
                << " does not contain patch " << thisPatch.name()
                << " in region " << pbm.mesh().name()
                << exit(FatalError);

            return -1;
        }


        if (patchIDs.size() == 2)
        {
            // Return the other patch
            return patchIDs[1-index];
        }
        else    // size == 1
        {
            return -1;
        }
    }
    else
    {
        if (patchIDs.size() != 1)
        {
            FatalErrorInFunction
                << "Couple patchGroup " << name()
                << " with contents " << patchIDs
                << " in region " << mesh.name()
                << " should only contain a single patch"
                << " when matching patch " << thisPatch.name()
                << " in region " << pbm.mesh().name()
                << exit(FatalError);
        }

        return patchIDs[0];
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::coupleGroupIdentifier::coupleGroupIdentifier()
:
    name_()
{}


Foam::coupleGroupIdentifier::coupleGroupIdentifier(const word& name)
:
    name_(name)
{}


Foam::coupleGroupIdentifier::coupleGroupIdentifier(const dictionary& dict)
:
    name_(dict.lookupOrDefault<word>("coupleGroup", ""))
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::label Foam::coupleGroupIdentifier::findOtherPatchID
(
    const polyPatch& thisPatch
) const
{
    const polyBoundaryMesh& pbm = thisPatch.boundaryMesh();

    return findOtherPatchID(pbm.mesh(), thisPatch);
}


Foam::label Foam::coupleGroupIdentifier::findOtherPatchID
(
    const polyPatch& thisPatch,
    word& otherRegion
) const
{
    const polyBoundaryMesh& pbm = thisPatch.boundaryMesh();
    const polyMesh& thisMesh = pbm.mesh();
    const Time& runTime = thisMesh.time();


    // Loop over all regions to find other patch in coupleGroup
    HashTable<const polyMesh*> meshSet = runTime.lookupClass<polyMesh>();

    label otherPatchID = -1;

    forAllConstIter(HashTable<const polyMesh*>, meshSet, iter)
    {
        const polyMesh& mesh = *iter();

        label patchID = findOtherPatchID(mesh, thisPatch);

        if (patchID != -1)
        {
            if (otherPatchID != -1)
            {
                FatalErrorInFunction
                    << "Couple patchGroup " << name()
                    << " should be present on only two patches"
                    << " in any of the meshes in " << meshSet.sortedToc()
                    << endl
                    << "    It seems to be present on patch "
                    << thisPatch.name()
                    << " in region " << thisMesh.name()
                    << ", on patch " << otherPatchID
                    << " in region " << otherRegion
                    << " and on patch " << patchID
                    << " in region " << mesh.name()
                    << exit(FatalError);
            }
            otherPatchID = patchID;
            otherRegion = mesh.name();
        }
    }

    if (otherPatchID == -1)
    {
        FatalErrorInFunction
            << "Couple patchGroup " << name()
            << " not found in any of the other meshes " << meshSet.sortedToc()
            << " on patch " << thisPatch.name()
            << " region " << thisMesh.name()
            << exit(FatalError);
    }

    return otherPatchID;
}


void Foam::coupleGroupIdentifier::write(Ostream& os) const
{
    if (valid())
    {
        writeEntry(os, "coupleGroup", name());
    }
}


// * * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const coupleGroupIdentifier& p)
{
    p.write(os);
    os.check("Ostream& operator<<(Ostream& os, const coupleGroupIdentifier& p");
    return os;
}


// ************************************************************************* //
