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

\*---------------------------------------------------------------------------*/

#include "stringOps.H"
#include "OSspecific.H"
#include "etcFiles.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Find the type/position of the ":-" or ":+" alternative values
static inline int findParameterAlternative
(
    const string& s,
    string::size_type& pos,
    string::size_type endPos
)
{
    while (pos != std::string::npos)
    {
        pos = s.find(':', pos);
        if (pos != std::string::npos)
        {
            if (pos < endPos)
            {
                // Nn-range: check for '+' or '-' following the ':'
                const int altType = s[pos+1];
                if (altType == '+' || altType == '-')
                {
                    return altType;
                }

                ++pos;    // Unknown/unsupported - continue at next position
            }
            else
            {
                // Out-of-range: abort
                pos = std::string::npos;
            }
        }
    }

    return 0;
}


//- Get dictionary or (optionally) environment variable
string getVariable
(
    const word& name,
    const dictionary& dict,
    const bool allowEnvVars,
    const bool allowEmpty
)
{
    const entry* ePtr = dict.lookupScopedEntryPtr(name, true, false);

    if (ePtr)
    {
        OStringStream buf;

        // Force floating point numbers to be printed with at least
        // some decimal digits.
        buf << scientific;

        buf.precision(IOstream::defaultPrecision());

        // Fail for non-primitiveEntry
        dynamicCast<const primitiveEntry>(*ePtr).write(buf, true);

        return buf.str();
    }
    else if (allowEnvVars)
    {
        string::const_iterator iter = name.begin();

        // Search for the sub-set of characters in name which are allowed
        // in environment variables
        string::size_type begVar = 0;
        string::size_type endVar = begVar;
        while (iter != name.end() && (isalnum(*iter) || *iter == '_'))
        {
            ++iter;
            ++endVar;
        }

        const word varName(name.substr(begVar, endVar - begVar), false);

        string varValue = getEnv(varName);

        if (!allowEmpty && varValue.empty())
        {
            FatalIOErrorInFunction
            (
                dict
            )   << "Cannot find dictionary or environment variable "
                << name << exit(FatalIOError);
        }

        varValue += name.substr(endVar, name.size() - endVar);

        return varValue;
    }
    else
    {
        FatalIOErrorInFunction
        (
            dict
        )   << "Cannot find dictionary variable "
            << name << exit(FatalIOError);

        return string::null;
    }
}


//- Recursively expands dictionary or environment variable starting at index
//  in string. Updates index.
string expand
(
    const string& s,
    string::size_type& index,
    const dictionary& dict,
    const bool allowEnvVars,
    const bool allowEmpty
)
{
    string newString;

    while (index < s.size())
    {
        if (s[index] == '$' && s[index+1] == '{')
        {
            // Recurse to parse variable name
            index += 2;
            string val = expand(s, index, dict, allowEnvVars, allowEmpty);
            newString.append(val);
        }
        else if (s[index] == '}')
        {
            return getVariable(newString, dict, allowEnvVars, allowEmpty);
        }
        else
        {
            newString.append(string(s[index]));
        }
        index++;
    }
    return newString;
}


//- Expand path parts of a string
Foam::string& inplaceExpandPath(string& s)
{
    if (!s.empty())
    {
        if (s[0] == '~')
        {
            // Expand initial ~
            //   ~/        => home directory
            //   ~OpenFOAM => site/user OpenFOAM configuration directory
            //   ~user     => home directory for specified user

            string user;
            fileName file;

            const string::size_type pos = s.find('/');
            if (pos != string::npos)
            {
                user = s.substr(1, pos - 1);
                file = s.substr(pos + 1);
            }
            else
            {
                user = s.substr(1);
            }

            // NB: be a bit lazy and expand ~unknownUser as an
            // empty string rather than leaving it untouched.
            // otherwise add extra test
            if (user == "OpenFOAM")
            {
                s = findEtcFile(file);
            }
            else
            {
                s = home(user)/file;
            }
        }
        else if (s[0] == '.')
        {
            // Expand a lone '.' and an initial './' into cwd
            if (s.size() == 1)
            {
                s = cwd();
            }
            else if (s[1] == '/')
            {
                s.std::string::replace(0, 1, cwd());
            }
        }
    }

    return s;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::string Foam::stringOps::expandEnvVar
(
    const string& original,
    const bool allowEmpty
)
{
    string s(original);
    return inplaceExpandEnvVar(s, allowEmpty);
}


Foam::string& Foam::stringOps::inplaceExpandEnvVar
(
    string& s,
    const bool allowEmpty
)
{
    string::size_type begVar = 0;

    // Expand $VARS
    // Repeat until nothing more is found
    while
    (
        (begVar = s.find('$', begVar)) != string::npos
     && begVar < s.size()-1
    )
    {
        if (begVar == 0 || s[begVar-1] != '\\')
        {
            // Find end of first occurrence
            string::size_type endVar = begVar;
            string::size_type delim = 0;

            // The type/position of the ":-" or ":+" alternative values
            int altType = 0;
            string::size_type altPos = string::npos;

            if (s[begVar+1] == '{')
            {
                endVar = s.find('}', begVar);
                delim = 1;

                // Check for ${parameter:-word} or ${parameter:+word}
                if (endVar != string::npos)
                {
                    altPos = begVar;
                    altType = findParameterAlternative(s, altPos, endVar);
                }
            }
            else
            {
                string::iterator iter = s.begin() + begVar + 1;

                while
                (
                    iter != s.end()
                 && (isalnum(*iter) || *iter == '_')
                )
                {
                    ++iter;
                    ++endVar;
                }
            }

            if (endVar == string::npos)
            {
                // Likely parsed '${...' without closing '}' - abort
                break;
            }
            else if (endVar == begVar)
            {
                // Parsed '${}' or $badChar  - skip over
                begVar = endVar + 1;
            }
            else
            {
                const word varName
                (
                    s.substr
                    (
                        begVar + 1 + delim,
                        (
                            (altPos == string::npos ? endVar : altPos)
                          - begVar - 2*delim
                        )
                    ),
                    false
                );

                std::string altValue;
                if (altPos != string::npos)
                {
                    // Had ":-" or ":+" alternative value
                    altValue = s.substr
                    (
                        altPos + 2,
                        endVar - altPos - 2*delim
                    );
                }

                const string varValue = getEnv(varName);
                if (varValue.size())
                {
                    if (altPos != string::npos && altType == '+')
                    {
                        // Was found, use ":+" alternative
                        s.std::string::replace
                        (
                            begVar,
                            endVar - begVar + 1,
                            altValue
                        );
                        begVar += altValue.size();
                    }
                    else
                    {
                        // Was found, use value
                        s.std::string::replace
                        (
                            begVar,
                            endVar - begVar + 1,
                            varValue
                        );
                        begVar += varValue.size();
                    }
                }
                else if (altPos != string::npos)
                {
                    // Use ":-" or ":+" alternative values
                    if (altType == '-')
                    {
                        // Was not found, use ":-" alternative
                        s.std::string::replace
                        (
                            begVar,
                            endVar - begVar + 1,
                            altValue
                        );
                        begVar += altValue.size();
                    }
                    else
                    {
                        // Was not found, ":+" alternative implies
                        // substitute with nothing
                        s.std::string::erase(begVar, endVar - begVar + 1);
                    }
                }
                else if (allowEmpty)
                {
                    s.std::string::erase(begVar, endVar - begVar + 1);
                }
                else
                {
                    FatalErrorInFunction
                        << "Unknown variable name '" << varName << "'"
                        << exit(FatalError);
                }
            }
        }
        else
        {
            ++begVar;
        }
    }

    return inplaceExpandPath(s);
}


Foam::string& Foam::stringOps::inplaceExpandCodeString
(
    string& s,
    const dictionary& dict,
    const word& dictVar,
    const char sigil
)
{
    string::size_type begVar = 0;

    // Expand $VAR or ${VAR}
    // Repeat until nothing more is found
    while
    (
        (begVar = s.find(sigil, begVar)) != string::npos
     && begVar < s.size()-1
    )
    {
        if (begVar == 0 || s[begVar-1] != '\\')
        {
            // Find end of first occurrence
            string::size_type endVar = begVar;
            string::size_type begDelim = 0, endDelim = 0;

            // Get the type, if any
            word varType;
            if (s[begVar+1] == '<')
            {
                begDelim += s.findClosing('>', begVar+1) - begVar;
                varType = s.substr(begVar+2, begDelim-2);
            }

            // Parse any braces and find the end of the variable
            if (s[begVar+begDelim+1] == '{')
            {
                endVar = s.find('}', begVar+begDelim);
                begDelim += 1;
                endDelim += 1;
            }
            else
            {
                endVar = begVar + begDelim;

                string::iterator iter = s.begin() + begVar + begDelim + 1;

                // Accept all dictionary and environment variable characters
                while
                (
                    iter != s.end()
                 &&
                    (
                        isalnum(*iter)
                     || *iter == '/' // For dictionary slash syntax
                     || *iter == '!' // For dictionary slash syntax
                     || *iter == '.' // For dictionary dot syntax
                     || *iter == ':' // For dictionary dot syntax
                     || *iter == '_'
                    )
                )
                {
                    ++iter;
                    ++endVar;
                }
            }

            if (endVar == string::npos)
            {
                // Likely parsed '${...' without closing '}' - abort
                break;
            }
            else if (endVar == begVar)
            {
                // Parsed '${}' or $badChar  - skip over
                begVar = endVar + 1;
            }
            else
            {
                const word varName
                (
                    s.substr
                    (
                        begVar + begDelim + 1,
                        (endVar - endDelim) - (begVar + begDelim)
                    ),
                    false
                );

                // Lookup in the dictionary
                const entry* ePtr = dict.lookupScopedEntryPtr
                (
                    varName,
                    true,
                    false   // Wildcards disabled. See primitiveEntry
                );

                // Substitute if found
                if (ePtr)
                {
                    OStringStream buf;

                    if (!dictVar.empty() && !varType.empty())
                    {
                        // If the dictionary is accessible and the type is
                        // known, then lookup the variable in the code, rather
                        // than substituting its value. That way we don't need
                        // to recompile this string if the value changes.
                        buf << dictVar
                            << ".lookupScoped<" << varType << ">"
                            << "(\"" << varName << "\", true, false)";
                    }

                    if (dictVar.empty() && !varType.empty())
                    {
                        // If the dictionary is not accessible but the type is
                        // known, then read the substituted value from a string
                        buf << "read<" << varType << ">(\"";
                    }

                    if (dictVar.empty() || varType.empty())
                    {
                        // If the dictionary is not accessible and/or the type
                        // is not known, then we need to substitute the
                        // variable's value

                        // Make sure floating point values print with at least
                        // some decimal points
                        buf << scientific;
                        buf.precision(IOstream::defaultPrecision());

                        // Write the dictionary or primitive entry. Fail if
                        // anything else.
                        if (ePtr->isDict())
                        {
                            ePtr->dict().write(buf, false);
                        }
                        else
                        {
                            dynamicCast<const primitiveEntry>(*ePtr)
                                .write(buf, true);
                        }
                    }

                    if (dictVar.empty() && !varType.empty())
                    {
                        // Close the string and read function as necessary
                        buf << "\")";
                    }

                    s.std::string::replace
                    (
                        begVar,
                        endVar - begVar + 1,
                        buf.str()
                    );
                    begVar += buf.str().size();
                }
                else
                {
                    // Not found. Leave original string untouched.
                    begVar = endVar + 1;
                }
            }
        }
        else
        {
            ++begVar;
        }
    }

    return s;
}


Foam::string& Foam::stringOps::inplaceExpandCodeTemplate
(
    string& s,
    const HashTable<string, word, string::hash>& mapping,
    const char sigil
)
{
    string::size_type begVar = 0;

    // Expand $VAR or ${VAR}
    // Repeat until nothing more is found
    while
    (
        (begVar = s.find(sigil, begVar)) != string::npos
     && begVar < s.size()-1
    )
    {
        if (begVar == 0 || s[begVar-1] != '\\')
        {
            // Find end of first occurrence
            string::size_type endVar = begVar;
            string::size_type delim = 0;

            if (s[begVar+1] == '{')
            {
                endVar = s.find('}', begVar);
                delim = 1;
            }
            else
            {
                string::iterator iter = s.begin() + begVar + 1;

                // Accept all dictionary and environment variable characters
                while
                (
                    iter != s.end()
                 &&
                    (
                        isalnum(*iter)
                     || *iter == '/' // For dictionary slash syntax
                     || *iter == '!' // For dictionary slash syntax
                     || *iter == '.' // For dictionary dot syntax
                     || *iter == ':' // For dictionary dot syntax
                     || *iter == '_'
                    )
                )
                {
                    ++iter;
                    ++endVar;
                }
            }

            if (endVar == string::npos)
            {
                // Likely parsed '${...' without closing '}' - abort
                break;
            }
            else if (endVar == begVar)
            {
                // Parsed '${}' or $badChar  - skip over
                begVar = endVar + 1;
            }
            else
            {
                const word varName
                (
                    s.substr
                    (
                        begVar + 1 + delim,
                        endVar - begVar - 2*delim
                    ),
                    false
                );

                // Lookup in the table
                HashTable<string, word, string::hash>::const_iterator fnd =
                    mapping.find(varName);

                // Substitute if found
                if (fnd != HashTable<string, word, string::hash>::end())
                {
                    s.std::string::replace
                    (
                        begVar,
                        endVar - begVar + 1,
                        *fnd
                    );
                    begVar += (*fnd).size();
                }
                else
                {
                    // Not found. Remove.
                    s.std::string::erase(begVar, endVar - begVar + 1);
                }
            }
        }
        else
        {
            ++begVar;
        }
    }

    return s;
}


Foam::string& Foam::stringOps::inplaceExpandEntry
(
    string& s,
    const dictionary& dict,
    const bool allowEnvVars,
    const bool allowEmpty,
    const char sigil
)
{
    string::size_type begVar = 0;

    // Expand $VAR or ${VAR}
    // Repeat until nothing more is found
    while
    (
        (begVar = s.find(sigil, begVar)) != string::npos
     && begVar < s.size()-1
    )
    {
        if (begVar == 0 || s[begVar-1] != '\\')
        {
            if (s[begVar+1] == '{')
            {
                // Recursive variable expansion mode
                label stringStart = begVar;
                begVar += 2;
                string varValue
                (
                    expand
                    (
                        s,
                        begVar,
                        dict,
                        allowEnvVars,
                        allowEmpty
                    )
                );

                s.std::string::replace
                (
                    stringStart,
                    begVar - stringStart + 1,
                    varValue
                );

                begVar = stringStart+varValue.size();
            }
            else
            {
                string::iterator iter = s.begin() + begVar + 1;

                // Accept all dictionary and environment variable characters
                string::size_type endVar = begVar;
                while
                (
                    iter != s.end()
                 &&
                    (
                        isalnum(*iter)
                     || *iter == '/' // For dictionary slash syntax
                     || *iter == '!' // For dictionary slash syntax
                     || *iter == '.' // For dictionary dot syntax
                     || *iter == ':' // For dictionary dot syntax
                     || *iter == '_'
                    )
                )
                {
                    ++iter;
                    ++endVar;
                }

                const word varName
                (
                    s.substr
                    (
                        begVar + 1,
                        endVar - begVar
                    ),
                    false
                );

                string varValue
                (
                    getVariable
                    (
                        varName,
                        dict,
                        allowEnvVars,
                        allowEmpty
                    )
                );

                s.std::string::replace
                (
                    begVar,
                    varName.size()+1,
                    varValue
                );
                begVar += varValue.size();
            }
        }
        else
        {
            ++begVar;
        }
    }

    return inplaceExpandPath(s);
}


Foam::string Foam::stringOps::trimLeft(const string& s)
{
    if (!s.empty())
    {
        string::size_type beg = 0;
        while (beg < s.size() && isspace(s[beg]))
        {
            ++beg;
        }

        if (beg)
        {
            return s.substr(beg);
        }
    }

    return s;
}


Foam::string& Foam::stringOps::inplaceTrimLeft(string& s)
{
    if (!s.empty())
    {
        string::size_type beg = 0;
        while (beg < s.size() && isspace(s[beg]))
        {
            ++beg;
        }

        if (beg)
        {
            s.erase(0, beg);
        }
    }

    return s;
}


Foam::string Foam::stringOps::trimRight(const string& s)
{
    if (!s.empty())
    {
        string::size_type sz = s.size();
        while (sz && isspace(s[sz-1]))
        {
            --sz;
        }

        if (sz < s.size())
        {
            return s.substr(0, sz);
        }
    }

    return s;
}


Foam::string& Foam::stringOps::inplaceTrimRight(string& s)
{
    if (!s.empty())
    {
        string::size_type sz = s.size();
        while (sz && isspace(s[sz-1]))
        {
            --sz;
        }

        s.resize(sz);
    }

    return s;
}


Foam::string Foam::stringOps::trim(const string& original)
{
    return trimLeft(trimRight(original));
}


Foam::string& Foam::stringOps::inplaceTrim(string& s)
{
    inplaceTrimRight(s);
    inplaceTrimLeft(s);

    return s;
}


Foam::string Foam::stringOps::breakIntoIndentedLines
(
    const string& message,
    const string::size_type nLength,
    const string::size_type nIndent
)
{
    const string indent(nIndent, token::SPACE);

    string result;

    word::size_type i0 = 0, i1 = 0;
    while (true)
    {
        const word::size_type iNewLine =
            message.find_first_of(token::NL, i1);
        const word::size_type iSpace =
            message.find_first_of(token::SPACE, i1);

        // New line next
        if
        (
            iNewLine != string::npos
         && (iSpace == string::npos || iNewLine < iSpace)
        )
        {
            result += indent + message.substr(i0, iNewLine - i0) + '\n';
            i0 = i1 = iNewLine + 1;
        }

        // Space next
        else if (iSpace != string::npos)
        {
            if (iSpace - i0 > nLength - nIndent)
            {
                result += indent + message.substr(i0, i1 - i0) + '\n';
                i0 = i1;
            }
            else
            {
                i1 = iSpace + 1;
            }
        }

        // End of string
        else
        {
            result += indent + message.substr(i0);
            break;
        }
    }

    return result;
}


// ************************************************************************* //
