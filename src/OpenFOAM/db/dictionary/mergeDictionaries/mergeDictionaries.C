/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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

#include "mergeDictionaries.H"
#include "ListOps.H"
#include "stringListOps.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

bool addEntry
(
    dictionary& tgtDict,
    entry& tgtEntry,
    const entry& srcEntry,
    const bool wildcards,
    const HashTable<wordList, word>& shortcuts
)
{
    bool changed = false;

    // Recursively merge sub-dictionaries
    // TODO: merge without copying
    if (tgtEntry.isDict() && srcEntry.isDict())
    {
        if
        (
            mergeDictionaries
            (
                const_cast<dictionary&>(tgtEntry.dict()),
                srcEntry.dict(),
                wildcards,
                shortcuts
            )
        )
        {
            changed = true;
        }
    }
    else
    {
        // Should use in-place modification instead of adding
        tgtDict.add(srcEntry.clone(tgtDict).ptr(), true);
        changed = true;
    }

    return changed;
}


labelList findMatches
(
    const bool wildcards,
    const HashTable<wordList, word>& shortcuts,
    const wordList& shortcutNames,
    const wordList& tgtKeys,
    const keyType& srcKey
)
{
    labelList matches;

    if (wildcards && srcKey.isPattern())
    {
        // Wildcard match
        matches = findStrings(srcKey, tgtKeys);
    }
    else if (shortcuts.size())
    {
        // See if shortcuts expand to valid tgtKeys
        labelList indices = findStrings(srcKey, shortcutNames);
        forAll(indices, i)
        {
            const word& name = shortcutNames[indices[i]];
            const wordList& keys = shortcuts[name];
            forAll(keys, j)
            {
                const label index = findIndex(tgtKeys, keys[j]);
                if (index != -1)
                {
                    matches.append(index);
                }
            }
        }
    }

    return matches;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

bool Foam::mergeDictionaries
(
    dictionary& tgtDict,
    const dictionary& srcDict,
    const bool wildcards,
    const HashTable<wordList, word>& shortcuts
)
{
    const wordList shortcutNames(shortcuts.toc());

    bool changed = false;

    // Save current (non-wildcard) keys before adding items.
    HashSet<word> tgtKeysSet;
    {
        List<keyType> keys = tgtDict.keys(false);
        forAll(keys, i)
        {
            tgtKeysSet.insert(keys[i]);
        }
    }

    // Pass 1. All literal matches
    forAllConstIter(IDLList<entry>, srcDict, srcIter)
    {
        const keyType& srcKey = srcIter().keyword();

        if (srcKey[0] == '~')
        {
            word eraseKey = srcKey(1, srcKey.size()-1);

            if (tgtDict.remove(eraseKey))
            {
                // Mark tgtDict entry as having been match for wildcard
                // handling later on.
                tgtKeysSet.erase(eraseKey);
            }

            changed = true;
        }
        else if (!wildcards || !(srcKey.isPattern() || shortcuts.found(srcKey)))
        {
            entry* entryPtr = tgtDict.lookupEntryPtr
            (
                srcKey,
                false,              // recursive
                false               // patternMatch
            );

            if (entryPtr)
            {
                // Mark tgtDict entry as having been match for wildcard
                // handling later on.
                tgtKeysSet.erase(entryPtr->keyword());

                if
                (
                    addEntry
                    (
                        tgtDict,
                       *entryPtr,
                        srcIter(),
                        wildcards,
                        shortcuts
                    )
                )
                {
                    changed = true;
                }
            }
            else
            {
                // not found - just add
                tgtDict.add(srcIter().clone(tgtDict).ptr());
                changed = true;
            }
        }
    }

    // Pass 2. Wildcard or shortcut matches (if any) on any non-match keys.
    if (tgtKeysSet.size() > 0)
    {
        // Pick up remaining dictionary entries
        wordList tgtKeys(tgtKeysSet.toc());

        forAllConstIter(IDLList<entry>, srcDict, srcIter)
        {
            const keyType& srcKey = srcIter().keyword();

            if (srcKey[0] == '~')
            {
                word eraseKey = srcKey(1, srcKey.size()-1);

                // List of indices into tgtKeys
                labelList matches
                (
                    findMatches
                    (
                        wildcards,
                        shortcuts,
                        shortcutNames,
                        tgtKeys,
                        eraseKey
                    )
                );

                // Remove all matches
                forAll(matches, i)
                {
                    tgtKeysSet.erase(tgtKeys[matches[i]]);
                }

                changed = true;
            }
            else
            {
                // List of indices into tgtKeys
                labelList matches
                (
                    findMatches
                    (
                        wildcards,
                        shortcuts,
                        shortcutNames,
                        tgtKeys,
                        srcKey
                    )
                );

                // Add all matches
                forAll(matches, i)
                {
                    const word& tgtKey = tgtKeys[matches[i]];

                    entry& tgtEntry = const_cast<entry&>
                    (
                        tgtDict.lookupEntry(tgtKey, false, false)
                    );

                    if
                    (
                        addEntry
                        (
                            tgtDict,
                            tgtEntry,
                            srcIter(),
                            wildcards,
                            HashTable<wordList, word>(0)
                        )
                    )
                    {
                        changed = true;
                    }
                }
            }
        }
    }

    return changed;
}


// ************************************************************************* //
