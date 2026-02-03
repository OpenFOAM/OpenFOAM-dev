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

Application
    foamUnits

Description
    Print information about dimensions and unit conversions

Usage
    \b foamUnits
    \b foamUnits \<unit\>
    \b foamUnits \<unit1\> \<unit2\>

    Options:
      - \par -all \n
        Print information for all available dimensions and units

Note
    This utility can be run with no arguments, one argument or two arguments.
    If no arguments are given this utility will print the names of all
    dimensions and units. If one argument is given then this is taken to be the
    name of a dimension or a unit and information will be printed regarding its
    relationship to the corresponding fundamental dimension or unit. If two
    arguments are given then these must be units (not dimensions) and
    conversions between them will be printed.

Example usage:
    foamUnits
    foamUnits -all
    foamUnits specificHeatCapacity
    foamUnits "(mol/cm^3)^-0.5/s" "(kmol/m^3)^-0.5/s"

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "unitConversion.H"
#include "stringOps.H"
#include "IOobject.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

string standardUnitName
(
    const wordList& dimensionUnitNames,
    const wordList& dimlessUnitNames,
    const unitConversion& unit
)
{
    string result;

    for (label i = 0; i < dimensionSet::nDimensions; ++ i)
    {
        const scalar e = unit[static_cast<dimensionSet::dimensionType>(i)];
        if (e == 0) continue;
        result.append(dimensionUnitNames[i]);
        if (e != 1) result.append("^" + name(e));
        result.append(" ");
    }
    for (label i = 0; i < unitConversion::nDimlessUnits; ++ i)
    {
        const scalar e = unit[static_cast<unitConversion::dimlessUnitType>(i)];
        if (e == 0) continue;
        result.append(dimlessUnitNames[i]);
        if (e != 1) result.append("^" + name(e));
        result.append(" ");
    }

    return result(result.size() - 1);
}


template<class Type>
bool stringIs(const string& str, const HashTable<Type>& types)
{
    forAllConstIter(typename HashTable<Type>, types, iter)
    {
        const size_t i0 = str.find(iter.key());
        const size_t i1 = i0 + iter.key().size();

        if (i0 == std::string::npos) continue;

        if (i0 > 0 && isalnum(str[i0 - 1])) continue;

        if (i1 < str.size() && isalnum(str[i1])) continue;

        return true;
    }

    return false;
}


bool isFundamental(const dimensionSet& dimension)
{
    label result = 0;

    for (label i = 0; i < dimensionSet::nDimensions; ++ i)
    {
        const dimensionSet::dimensionType t =
            static_cast<dimensionSet::dimensionType>(i);

        result =
            result == 0 && dimension[t] == 1 ? 1
          : dimension[t] == 0 ? result
          : -1;
    }

    return result == 1;
}


bool isFundamental(const unitConversion& unit)
{
    label result = 0;

    for (label i = 0; i < dimensionSet::nDimensions; ++ i)
    {
        const dimensionSet::dimensionType t =
            static_cast<dimensionSet::dimensionType>(i);

        result =
            result == 0 && unit[t] == 1 ? 1
          : unit[t] == 0 ? result
          : -1;
    }

    for (label i = 0; i < unitConversion::nDimlessUnits; ++ i)
    {
        const unitConversion::dimlessUnitType t =
            static_cast<unitConversion::dimlessUnitType>(i);

        result =
            result == 0 && unit[t] == 1 ? 1
          : unit[t] == 0 ? result
          : -1;
    }

    return result == 1;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "removeCaseOptions.H"
    writeInfoHeader = false;

    argList::addBoolOption
    (
        "all",
        "list only the names of the dimensions and units"
    );

    const label nArgs = argList::nArgs(argc, argv);
    switch (nArgs)
    {
        case 1:
            argList::validArgs.append("unit");
            break;
        case 2:
            argList::validArgs.append("unit1");
            argList::validArgs.append("unit2");
            break;
        default:
            break;
    }

    argList::addNote
    (
        "This utility can be run with no arguments, one argument or two\n"
        "arguments. If no arguments are given this utility will print the\n"
        "names of all dimensions and units. If one argument is given then \n"
        "this is taken to be the name of a dimension or a unit and\n"
        "information will be printed regarding its relationship to the\n"
        "corresponding fundamental dimension or unit. If two arguments are\n"
        "given then these must be units (not dimensions) and conversions\n"
        "between them will be printed.\n"
        "\n"
        "Example usage:\n"
        "    foamUnits\n"
        "    foamUnits -all\n"
        "    foamUnits specificHeatCapacity\n"
        "    foamUnits \"(mol/cm^3)^-0.5/s\" \"(kmol/m^3)^-0.5/s\""
    );

    argList args(argc, argv);

    // Build lists of fundamental unit names
    wordList dimensionUnitNames(dimensionSet::nDimensions);
    wordList dimlessUnitNames(unitConversion::nDimlessUnits);
    forAllConstIter(HashTable<unitConversion>, units(), iter)
    {
        const unitConversion& unit = iter();

        label dimensioni = -1, dimlessUniti = -1;
        for (label i = 0; i < dimensionSet::nDimensions; ++ i)
        {
            const dimensionSet::dimensionType t =
                static_cast<dimensionSet::dimensionType>(i);
            if (dimensioni >= 0 && unit[t] != 0)
            {
                dimensioni = -2;
            }
            if (dimensioni == -1 && unit[t] == 1 && unit.standard())
            {
                dimensioni = i;
            }
        }
        for (label i = 0; i < unitConversion::nDimlessUnits; ++ i)
        {
            const unitConversion::dimlessUnitType t =
                static_cast<unitConversion::dimlessUnitType>(i);
            if (dimlessUniti >= 0 && unit[t] != 0)
            {
                dimlessUniti = -1;
            }
            if (dimlessUniti == -1 &&  unit[t] == 1 && unit.standard())
            {
                dimlessUniti = i;
            }
        }

        if (dimensioni >= 0 && dimlessUniti == -1)
        {
            dimensionUnitNames[dimensioni] = iter.key();
        }
        if (dimensioni == -1 && dimlessUniti >= 0)
        {
            dimlessUnitNames[dimlessUniti] = iter.key();
        }
    }

    auto printDimension = [&](const word& name, const dimensionSet& dimension)
    {
        Info<< "Dimension [" << name << "]" << nl;
        if (!isFundamental(dimension))
            Info<< "+ Dimensions = " << dimension.info() << nl;
        Info << "+ Exponents = " << dimension << nl
             << endl;
    };

    auto printUnit = [&](const word& name, const unitConversion& unit)
    {
        const string standardName =
            standardUnitName(dimensionUnitNames, dimlessUnitNames, unit);
        Info<< "Unit [" << name << "]" << nl
            << "+ Dimensions = " << unit.dimensions().info() << nl;
        if (!isFundamental(unit))
            Info<< "+ Standard Unit = [" << standardName.c_str() << "]" << nl;
        Info<< "+ Conversion Factor = " << unit.toStandard(scalar(1)) << nl
            << endl;
    };

    Info<< endl;

    // Print all dimensions and units
    if (args.optionFound("all"))
    {
        forAllConstIter(HashTable<dimensionSet>, dimensions(), iter)
        {
            printDimension(iter.key(), iter());
        }

        forAllConstIter(HashTable<unitConversion>, units(), iter)
        {
            printUnit(iter.key(), iter());
        }

        if (nArgs)
        {
            IOobject::writeDivider(Info);
            Info<< nl;
        }
    }

    // Print a single dimension or unit
    if (nArgs == 1)
    {
        const string name(args[1]);

        const bool isDimension = stringIs(name, dimensions());
        const bool isUnit = stringIs(name, units());

        if (isDimension && !isUnit)
        {
            printDimension(name, IStringStream(("[" + name + "]").c_str())());

            return 0;
        }

        if (!isDimension && isUnit)
        {
            printUnit(name, IStringStream(("[" + name + "]").c_str())());

            return 0;
        }

        FatalErrorInFunction
            << "'" << args[1].c_str()
            << "' is not a valid dimension or unit"
            << exit(FatalError);
    }

    // Print the conversion between two units
    if (nArgs == 2)
    {
        const string name1(args[1]);
        const string name2(args[2]);

        auto assertStringIsUnit = [](const string& str)
        {
            if (stringIs(str, dimensions()))
            {
                FatalErrorInFunction
                    << "'" << str.c_str() << "' is a dimension. "
                    << "Comparison is only supported for units."
                    << exit(FatalError);
            }
        };
        assertStringIsUnit(name1);
        assertStringIsUnit(name2);

        const unitConversion unit1
        (
            IStringStream(("[" + name1 + "]").c_str())()
        );
        const unitConversion unit2
        (
            IStringStream(("[" + name2 + "]").c_str())()
        );

        // Check the units are the same, except for the multiplier
        unitConversion
        (
            unit1.dimensions(),
            unit1[unitConversion::FRACTION],
            unit1[unitConversion::ANGLE],
            1
        )
      + unitConversion
        (
            unit2.dimensions(),
            unit2[unitConversion::FRACTION],
            unit2[unitConversion::ANGLE],
            1
        );

        const string standardName =
            standardUnitName(dimensionUnitNames, dimlessUnitNames, unit1);

        Info<< "Units [" << name1.c_str() << "] [" << name2.c_str() << ']' << nl
            << "+ Dimensions = " << unit1.dimensions().info() << nl
            << "+ Standard Unit = [" << standardName.c_str() << "]" << nl
            << "+ Conversion: " << 1 << " [" << name1.c_str() << "] = "
            << unit2.toUser(unit1.toStandard(scalar(1))) << " ["
            << name2.c_str() << ']' << nl
            << "+ Conversion: " << 1 << " [" << name2.c_str() << "] = "
            << unit1.toUser(unit2.toStandard(scalar(1))) << " ["
            << name1.c_str() << ']' << nl
            << endl;

        return 0;
    }

    // Print just the names of the dimensions and units
    if (!args.optionFound("all"))
    {
        auto print = [](const char* group, const DynamicList<word>& names)
        {
            Info<< group << ":" << nl;

            OStringStream str;
            forAll(names, i)
            {
                str << (i ? " " : "") << "[" << names[i] << ']';
            }

            Info<< stringOps::breakIntoIndentedLines(str.str(), 80, 4).c_str()
                << nl << endl;
        };

        DynamicList<word> fundamentalDimensions;
        DynamicList<word> derivedDimensions;
        forAllConstIter(HashTable<dimensionSet>, dimensions(), iter)
        {
            (
                isFundamental(iter())
              ? fundamentalDimensions
              : derivedDimensions
            ).append(iter.key());
        }

        print("Fundamental Dimensions", fundamentalDimensions);
        print("Derived Dimensions", derivedDimensions);

        DynamicList<word> fundamentalUnits;
        DynamicList<word> derivedUnits;
        DynamicList<word> scaledUnits;
        DynamicList<word> derivedScaledUnits;
        forAllConstIter(HashTable<unitConversion>, units(), iter)
        {
            (
                isFundamental(iter()) && iter().standard() ? fundamentalUnits
              : !isFundamental(iter()) && iter().standard() ? derivedUnits
              : isFundamental(iter()) && !iter().standard() ? scaledUnits
              : derivedScaledUnits
            ).append(iter.key());
        }

        print("Fundamental Units", fundamentalUnits);
        print("Derived Units", derivedUnits);
        print("Scaled Units", scaledUnits);
        print("Derived-Scaled Units", derivedScaledUnits);
    }

    return 0;
}

// ************************************************************************* //
