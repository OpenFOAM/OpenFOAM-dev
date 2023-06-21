/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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

Namespace
    Foam::stringOps

Description
    Collection of static functions to do various simple string-related
    operations

SourceFiles
    stringOps.C

\*---------------------------------------------------------------------------*/

#ifndef stringOps_H
#define stringOps_H

#include "string.H"
#include "dictionary.H"
#include "HashTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Namespace stringOps Declaration
\*---------------------------------------------------------------------------*/

namespace stringOps
{
    //- Expand all occurrences of environment variables and paths
    //  Expansion includes:
    //  -# environment variables
    //    - "$VAR", "${VAR}"
    //  -# current directory
    //    - leading "./" : the current directory
    //  -# tilde expansion
    //    - leading "~/" : home directory
    //    - leading "~user" : home directory for specified user
    //    - leading "~OpenFOAM" : site/user OpenFOAM configuration directory
    //
    //  Supports default values as per the Bourne/Korn shell.
    //  \code
    //      "${parameter:-defValue}"
    //  \endcode
    //  If parameter is unset or null, the \c defValue is substituted.
    //  Otherwise, the value of parameter is substituted.
    //
    //  Supports alternative values as per the Bourne/Korn shell.
    //  \code
    //      "${parameter:+altValue}"
    //  \endcode
    //  If parameter is unset or null, nothing is substituted.
    //  Otherwise the \c altValue is substituted.
    //
    //  Any unknown entries are removed silently, if allowEmpty is true.
    //
    //  Malformed entries (eg, brace mismatch, sigil followed by bad character)
    //  are left as is.
    //
    //  \sa
    //  Foam::findEtcFile
    string expandEnvVar
    (
        const string&,
        const bool allowEmpty = false
    );

    //- Expand all occurrences of environment variables and paths
    //  See expandEnvVar
    string& inplaceExpandEnvVar
    (
        string&,
        const bool allowEmpty = false
    );

    //- Inplace expand occurrences of variables according to the dictionary
    //  Expansion includes:
    //  -# variables
    //    - "$VAR", "${VAR}"
    //
    //  Any unknown entries are left as-is
    //
    //  \note the leading sigil can be changed to avoid conflicts with other
    //  string expansions
    string& inplaceExpandCodeString
    (
        string&,
        const dictionary& dict,
        const word& dictVar = "dict",
        const char sigil = '$'
    );

    //- Inplace expand occurrences of variables according to the mapping
    //  Expansion includes:
    //  -# variables
    //    - "$VAR", "${VAR}"
    //
    //  Any unknown entries are removed silently.
    //
    //  \note the leading sigil can be changed to avoid conflicts with other
    //  string expansions
    string& inplaceExpandCodeTemplate
    (
        string&,
        const HashTable<string, word, string::hash>& mapping,
        const char sigil = '$'
    );

    //- Inplace expand occurrences of variables according to the dictionary
    //  and optionally environment variables
    //  Expansion includes:
    //  -# variables
    //    - "$VAR", "${VAR}"
    //  -# current directory
    //    - leading "./" : the current directory
    //  -# tilde expansion
    //    - leading "~/" : home directory
    //    - leading "~user" : home directory for specified user
    //    - leading "~OpenFOAM" : site/user OpenFOAM configuration directory
    //
    //  with the "${}" syntax doing a recursive substitution.
    //
    //  Any unknown entries are left as-is
    //
    //  \note the leading sigil can be changed to avoid conflicts with other
    //  string expansions
    string& inplaceExpandEntry
    (
        string& s,
        const dictionary& dict,
        const bool allowEnvVars,
        const bool allowEmpty,
        const char sigil = '$'
    );


    //- Return string trimmed of leading whitespace
    string trimLeft(const string&);

    //- Trim leading whitespace inplace
    string& inplaceTrimLeft(string&);

    //- Return string trimmed of trailing whitespace
    string trimRight(const string&);

    //- Trim trailing whitespace inplace
    string& inplaceTrimRight(string&);

    //- Return string trimmed of leading and trailing whitespace
    string trim(const string&);

    //- Trim leading and trailing whitespace inplace
    string& inplaceTrim(string&);


    //- Break a string up into indented lines
    string breakIntoIndentedLines
    (
        const string& str,
        const string::size_type nLength = 80,
        const string::size_type nIndent = 0
    );


} // End namespace stringOps


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
