/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2022 OpenFOAM Foundation
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

#include "phaseInterface.H"
#include "phaseSystem.H"
#include "surfaceTensionModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    wordList phaseInterface::headSeparators_ = Foam::wordList();

    bool phaseInterfaceAddedHeadSeparator =
        phaseInterface::addHeadSeparator(Foam::word::null);

    HashTable<word> phaseInterface::oldSeparatorToSeparator_ =
        Foam::HashTable<word>();
}

namespace Foam
{
    defineTypeNameAndDebugWithName
    (
       phaseInterface,
       separatorsToTypeName(wordList(1, word::null)).c_str(),
       0
    );
    defineRunTimeSelectionTable(phaseInterface, word);
    addToRunTimeSelectionTable(phaseInterface, phaseInterface, word);
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

const Foam::phaseModel& Foam::phaseInterface::getPhase1
(
    const phaseModel& phase1,
    const phaseModel& phase2
)
{
    return phase1.index() < phase2.index() ? phase1 : phase2;
}


const Foam::phaseModel& Foam::phaseInterface::getPhase2
(
    const phaseModel& phase1,
    const phaseModel& phase2
)
{
    return phase1.index() < phase2.index() ? phase2 : phase1;
}


bool Foam::phaseInterface::addHeadSeparator(const word& separator)
{
    if (findIndex(headSeparators_, separator) == -1)
    {
        headSeparators_.append(separator);
        return true;
    }
    else
    {
        return false;
    }
}


bool Foam::phaseInterface::addOldSeparatorToSeparator
(
    const word& oldSeparator,
    const word& separator
)
{
    return oldSeparatorToSeparator_.insert(oldSeparator, separator);
}


Foam::wordList Foam::phaseInterface::nameToNameParts
(
    const phaseSystem& fluid,
    const word& name
)
{
    auto error = [&]()
    {
        FatalErrorInFunction
            << "Could not parse interface name \""
            << name << "\"" << exit(FatalError);
    };

    wordList nameParts;
    word::size_type i = 0;
    while (true)
    {
        // Match potentially multiple phases
        label nPhaseNameParts = 0;
        while (true)
        {
            // Match the longest possible phase name
            word phaseNamePart;
            forAll(fluid.phases(), phasei)
            {
                const word& phaseName =
                    fluid.phases()[phasei].name();
                if
                (
                    phaseNamePart.size() < phaseName.size()
                 && name.size() - i >= phaseName.size()
                 && name(i, phaseName.size()) == phaseName
                )
                {
                    nPhaseNameParts ++;
                    phaseNamePart = phaseName;
                }
            }
            if (phaseNamePart == word::null) break;

            // Add the phase name part
            nameParts.append(phaseNamePart);
            i += phaseNamePart.size();

            // Break if at the end
            if (i == name.size()) break;

            // Add an empty separator
            nameParts.append(word::null);

            // Pass the underscore
            if (name[i] != '_') error();
            i += 1;
        }

        // Error if we haven't matched any phases
        if (nPhaseNameParts == 0) error();

        // Break if at the end
        if (i == name.size()) break;

        // Match and add a non-empty separator
        const word::size_type j = name.find_first_of('_', i);
        if (j == word::npos) error();
        nameParts.last() = name(i, j - i);
        i = j;

        // Pass the underscore
        if (name[i] != '_') error();
        i += 1;
    }

    return nameParts;
}


Foam::wordList Foam::phaseInterface::nameToSeparators
(
    const phaseSystem& fluid,
    const word& name
)
{
    const wordList nameParts = nameToNameParts(fluid, name);

    wordList separators(nameParts.size()/2);
    forAll(separators, separatori)
    {
        separators[separatori] = nameParts[2*separatori + 1];
    }

    return separators;
}


Foam::word Foam::phaseInterface::separatorsToTypeName
(
    const wordList& separatorsUnsorted
)
{
    // Take a copy of the separators
    wordList separators(separatorsUnsorted);

    // Sort all but the first
    SubList<word> tailSeparators(separators, separators.size() - 1, 1);
    Foam::sort(tailSeparators);

    // Stitch back together using a placeholder phase name
    static const word phaseName = "<phase>";
    word typeName = phaseName;
    forAll(separators, separatori)
    {
        typeName.append
        (
            '_'
          + separators[separatori]
          + (separators[separatori].empty() ? "" : "_")
          + phaseName
        );
    }

    return typeName;
}


Foam::word Foam::phaseInterface::nameToTypeName
(
    const phaseSystem& fluid,
    const word& name
)
{
    return separatorsToTypeName(nameToSeparators(fluid, name));
}


Foam::word Foam::phaseInterface::namePartsToName
(
    const phaseSystem& fluid,
    const wordList& nameParts
)
{
    word name;

    forAll(nameParts, i)
    {
        if (nameParts[i] != word::null)
        {
            name.append(nameParts[i] + "_");
        }
    }

    return name(name.size() - 1);
}


Foam::word Foam::phaseInterface::oldNamePartsToName
(
    const phaseSystem& fluid,
    const wordList& oldNameParts
)
{
    word name;

    forAll(oldNameParts, i)
    {
        // Phase name part
        if (fluid.phases().found(oldNameParts[i]))
        {
            name.append(oldNameParts[i] + "_");
        }

        // Separator
        else
        {
            const word& oldSeparator = oldNameParts[i];
            const word& separator = oldSeparatorToSeparator_[oldSeparator];

            if (separator != word::null)
            {
                name.append(separator + "_");
            }
        }
    }

    return name(name.size() - 1);
}


Foam::Tuple2<const Foam::phaseModel&, const Foam::phaseModel&>
Foam::phaseInterface::identifyPhases
(
    const phaseSystem& fluid,
    const word& name,
    const wordList& separators
)
{
    const wordList nameParts = nameToNameParts(fluid, name);

    label nameParti = -1;
    bool multiple = false;
    for (label namePartj = 1; namePartj < nameParts.size() - 1; namePartj += 2)
    {
        forAll(separators, separatori)
        {
            if (nameParts[namePartj] == separators[separatori])
            {
                multiple = multiple || nameParti != -1;
                nameParti = namePartj;
            }
        }
    }

    if (nameParti == -1)
    {
        FatalErrorInFunction
            << "No matches identified in \"" << name
            << "\" for separators " << separators << exit(FatalError);
    }

    if (multiple)
    {
        FatalErrorInFunction
            << "Multiple matches identified in \"" << name
            << "\" for separators " << separators << exit(FatalError);
    }

    return
        Tuple2<const phaseModel&, const phaseModel&>
        (
            fluid.phases()[nameParts[nameParti - 1]],
            fluid.phases()[nameParts[nameParti + 1]]
        );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseInterface::phaseInterface
(
    const phaseModel& phase1,
    const phaseModel& phase2
)
:
    phase1_(getPhase1(phase1, phase2)),
    phase2_(getPhase2(phase1, phase2)),
    g_(phase1.mesh().lookupObject<uniformDimensionedVectorField>("g"))
{}


Foam::phaseInterface::phaseInterface
(
    const Tuple2<const phaseModel&, const phaseModel&>& phases
)
:
    phaseInterface(phases.first(), phases.second())
{}


Foam::phaseInterface::phaseInterface
(
    const phaseSystem& fluid,
    const word& name
)
:
    phaseInterface(identifyPhases(fluid, name, headSeparators_))
{}


Foam::phaseInterface::phaseInterface
(
    const phaseSystem& fluid,
    const phaseInterfaceKey& key
)
:
    phaseInterface(fluid.phases()[key.first()], fluid.phases()[key.second()])
{}


Foam::autoPtr<Foam::phaseInterface> Foam::phaseInterface::clone() const
{
    return New(fluid(), name());
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::phaseInterface::~phaseInterface()
{}


// * * * * * * * * * * * * * * * * Selector  * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::phaseInterface> Foam::phaseInterface::New
(
    const phaseSystem& fluid,
    const word& name
)
{
    wordConstructorTable::iterator cstrIter =
        wordConstructorTablePtr_->find(nameToTypeName(fluid, name));

    if (cstrIter == wordConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown phaseInterface type "
            << name << endl << endl
            << "Valid phaseInterface types are : " << endl
            << wordConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return cstrIter()(fluid, name);
}


Foam::autoPtr<Foam::phaseInterface> Foam::phaseInterface::New
(
    const phaseInterface& interface1,
    const phaseInterface& interface2
)
{
    const phaseSystem& fluid = interface1.fluid();

    auto error = [&]()
    {
        FatalErrorInFunction
            << "Could not combine interfaces " << interface1.name()
            << " and " << interface2.name() << exit(FatalError);
    };

    // Split the names into parts
    const wordList nameParts1 = nameToNameParts(fluid, interface1.name());
    const wordList nameParts2 = nameToNameParts(fluid, interface2.name());

    // Check the heads are consistent
    const Pair<word> headNames1(nameParts1[0], nameParts1[2]);
    const Pair<word> headNames2(nameParts2[0], nameParts2[2]);
    const word& headSeparator1 = nameParts1[1];
    const word& headSeparator2 = nameParts2[1];
    const label headNamesCompare = Pair<word>::compare(headNames1, headNames2);
    const bool haveBothHeadSeparators =
        headSeparator1 != word::null && headSeparator2 != word::null;
    if
    (
        (headNamesCompare == 0)
     || (haveBothHeadSeparators && headSeparator1 != headSeparator2)
     || (haveBothHeadSeparators && headNamesCompare != 1)
    )
    {
        error();
    }

    // Add the head to the list
    wordList nameParts
    (
        SubList<word>
        (
            headSeparator1 != word::null ? nameParts1 : nameParts2,
            3
        )
    );

    // Add the tail from interface 1
    for (label i1 = 3; i1 < nameParts1.size(); i1 += 2)
    {
        nameParts.append(nameParts1[i1]);
        nameParts.append(nameParts1[i1 + 1]);
    }

    // Add the tail from interface 2, filtering out duplicates
    for (label i2 = 3; i2 < nameParts2.size(); i2 += 2)
    {
        bool append2 = true;
        for (label i1 = 3; i1 < nameParts1.size(); i1 += 2)
        {
            if (nameParts1[i1] == nameParts2[i2])
            {
                if (nameParts1[i1 + 1] != nameParts2[i2 + 1]) error();
                append2 = false;
            }
        }
        if (append2)
        {
            nameParts.append(nameParts2[i2]);
            nameParts.append(nameParts2[i2 + 1]);
        }
    }

    // Select the new combined interface
    return New(fluid, namePartsToName(fluid, nameParts));
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::word Foam::phaseInterface::name() const
{
    return phase1().name() + "_" + phase2().name();
}


Foam::tmp<Foam::volScalarField> Foam::phaseInterface::rho() const
{
    return phase1()*phase1().rho() + phase2()*phase2().rho();
}


Foam::tmp<Foam::volScalarField> Foam::phaseInterface::magUr() const
{
    return mag(phase1().U() - phase2().U());
}


Foam::tmp<Foam::volScalarField> Foam::phaseInterface::sigma() const
{
    return fluid().sigma(*this);
}


// ************************************************************************* //
