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

#include "polyTopoChanger.H"
#include "polyMesh.H"
#include "polyTopoChange.H"
#include "Time.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(polyTopoChanger, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

void Foam::polyTopoChanger::readModifiers()
{
    if
    (
        readOpt() == IOobject::MUST_READ
     || readOpt() == IOobject::MUST_READ_IF_MODIFIED
     || (readOpt() == IOobject::READ_IF_PRESENT && headerOk())
    )
    {
        if (readOpt() == IOobject::MUST_READ_IF_MODIFIED)
        {
            WarningInFunction
                << "Specified IOobject::MUST_READ_IF_MODIFIED but class"
                << " does not support automatic re-reading."
                << endl;
        }


        PtrList<polyMeshModifier>& modifiers = *this;

        // Read modifiers
        Istream& is = readStream(typeName);

        PtrList<entry> patchEntries(is);
        modifiers.setSize(patchEntries.size());

        forAll(modifiers, modifierI)
        {
            modifiers.set
            (
                modifierI,
                polyMeshModifier::New
                (
                    patchEntries[modifierI].keyword(),
                    patchEntries[modifierI].dict(),
                    modifierI,
                    *this
                )
            );
        }

        // Check state of IOstream
        is.check("polyTopoChanger::readModifiers()");

        close();
    }
}


Foam::polyTopoChanger::polyTopoChanger
(
    const IOobject& io,
    polyMesh& mesh
)
:
    PtrList<polyMeshModifier>(),
    regIOobject(io),
    mesh_(mesh)
{
    readModifiers();
}


Foam::polyTopoChanger::polyTopoChanger(polyMesh& mesh)
:
    PtrList<polyMeshModifier>(),
    regIOobject
    (
        IOobject
        (
            "meshModifiers",
            mesh.time().findInstance
            (
                mesh.meshDir(),
                "meshModifiers",
                IOobject::READ_IF_PRESENT
            ),
            mesh.meshSubDir,
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        )
    ),
    mesh_(mesh)
{
    readModifiers();
}


Foam::wordList Foam::polyTopoChanger::types() const
{
    const PtrList<polyMeshModifier>& modifiers = *this;

    wordList t(modifiers.size());

    forAll(modifiers, modifierI)
    {
        t[modifierI] = modifiers[modifierI].type();
    }

    return t;
}


Foam::wordList Foam::polyTopoChanger::names() const
{
    const PtrList<polyMeshModifier>& modifiers = *this;

    wordList t(modifiers.size());

    forAll(modifiers, modifierI)
    {
        t[modifierI] = modifiers[modifierI].name();
    }

    return t;
}


bool Foam::polyTopoChanger::changeTopology() const
{
    // Go through all mesh modifiers and accumulate the morphing information
    const PtrList<polyMeshModifier>& topoChanges = *this;

    bool triggerChange = false;

    forAll(topoChanges, morphI)
    {
        if (topoChanges[morphI].active())
        {
            bool curTriggerChange = topoChanges[morphI].changeTopology();

            if (debug)
            {
                Info<< "Modifier " << morphI << " named "
                    << topoChanges[morphI].name();

                if (curTriggerChange)
                {
                    Info<< " morphing" << endl;
                }
                else
                {
                    Info<< " unchanged" << endl;
                }
            }

            triggerChange = triggerChange || curTriggerChange;
        }
        else
        {
            if (debug)
            {
                Info<< "Modifier " << morphI  << " named "
                    << topoChanges[morphI].name() << " inactive" << endl;
            }
        }

    }

    return triggerChange;
}


Foam::autoPtr<Foam::polyTopoChange>
Foam::polyTopoChanger::topoChangeRequest() const
{
    // Collect changes from all modifiers
    const PtrList<polyMeshModifier>& topoChanges = *this;

    polyTopoChange* refPtr(new polyTopoChange(mesh()));
    polyTopoChange& ref = *refPtr;

    forAll(topoChanges, morphI)
    {
        if (topoChanges[morphI].active())
        {
            topoChanges[morphI].setRefinement(ref);
        }
    }

    return autoPtr<polyTopoChange>(refPtr);
}


void Foam::polyTopoChanger::modifyMotionPoints(pointField& p) const
{
    const PtrList<polyMeshModifier>& topoChanges = *this;

    forAll(topoChanges, morphI)
    {
        if (topoChanges[morphI].active())
        {
            topoChanges[morphI].modifyMotionPoints(p);
        }
    }
}


void Foam::polyTopoChanger::update(const mapPolyMesh& m)
{
    // Go through all mesh modifiers and accumulate the morphing information
    PtrList<polyMeshModifier>& topoChanges = *this;

    forAll(topoChanges, morphI)
    {
        topoChanges[morphI].updateMesh(m);
    }

    // Force the mesh modifiers to auto-write.  This allows us to
    // preserve the current state of modifiers corresponding with
    // the mesh.
    writeOpt() = IOobject::AUTO_WRITE;
    instance() = mesh_.time().timeName();
}


Foam::autoPtr<Foam::mapPolyMesh> Foam::polyTopoChanger::changeMesh
(
    const bool inflate,
    const bool syncParallel,
    const bool orderCells,
    const bool orderPoints
)
{
    if (changeTopology())
    {
        autoPtr<polyTopoChange> ref = topoChangeRequest();

        autoPtr<mapPolyMesh> topoChangeMap = ref().changeMesh
        (
            mesh_,
            inflate,
            syncParallel,
            orderCells,
            orderPoints
        );

        update(topoChangeMap());
        mesh_.updateMesh(topoChangeMap());
        return topoChangeMap;
    }
    else
    {
        return autoPtr<mapPolyMesh>(nullptr);
    }
}


void Foam::polyTopoChanger::addTopologyModifiers
(
    const List<polyMeshModifier*>& tm
)
{
    setSize(tm.size());

    // Copy the patch pointers
    forAll(tm, tmI)
    {
        if (tm[tmI]->topoChanger() != *this)
        {
            FatalErrorInFunction
                << "Mesh modifier created with different mesh reference."
                << abort(FatalError);
        }
        set(tmI, tm[tmI]);
    }

    writeOpt() = IOobject::AUTO_WRITE;
}


Foam::label Foam::polyTopoChanger::findModifierID
(
    const word& modName
) const
{
    const PtrList<polyMeshModifier>& topoChanges = *this;

    forAll(topoChanges, morphI)
    {
        if (topoChanges[morphI].name() == modName)
        {
            return morphI;
        }
    }

    // Modifier not found
    if (debug)
    {
        WarningInFunction
            << "List of available modifier names: " << names() << endl;
    }

    // Not found, return -1
    return -1;
}


bool Foam::polyTopoChanger::writeData(Ostream& os) const
{
    os << *this;
    return os.good();
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

bool Foam::polyTopoChanger::operator!=(const polyTopoChanger& me) const
{
    return &me != this;
}


bool Foam::polyTopoChanger::operator==(const polyTopoChanger& me) const
{
    return &me == this;
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const polyTopoChanger& mme)
{
    os  << mme.size() << nl << token::BEGIN_LIST;

    forAll(mme, mmeI)
    {
        mme[mmeI].writeDict(os);
    }

    os  << token::END_LIST;

    return os;
}


// ************************************************************************* //
