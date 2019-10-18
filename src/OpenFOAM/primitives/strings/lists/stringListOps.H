/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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

InNamespace
    Foam

Description
    Operations on lists of strings.

SourceFiles
    stringListOpsTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef stringListOps_H
#define stringListOps_H

#include "regExp.H"
#include "labelList.H"
#include "stringList.H"
#include "wordReList.H"
#include "wordReListMatcher.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    // single-string matches:

    //- Return true if string matches one of the regular expressions
    inline bool findStrings
    (
        const wordReListMatcher& matcher,
        const std::string& str
    )
    {
        return matcher.match(str);
    }

    // multi-string matches:

    //- Return list indices for matching strings
    template<class Matcher, class StringType>
    labelList findMatchingStrings
    (
        const Matcher&,
        const UList<StringType>&,
        const bool invert=false
    );

    //- Return list indices for strings matching the regular expression
    //  Template partial specialization of findMatchingStrings
    template<class StringType>
    labelList findStrings
    (
        const regExp& re,
        const UList<StringType>& lst,
        const bool invert=false
    )
    {
        return findMatchingStrings(re, lst, invert);
    }

    //- Return list indices for strings matching the regular expression
    //  Template partial specialization of findMatchingStrings
    template<class StringType>
    labelList findStrings
    (
        const char* rePattern,
        const UList<StringType>& lst,
        const bool invert=false
    )
    {
        const regExp re(rePattern);
        return findStrings(re, lst, invert);
    }

    //- Return list indices for strings matching the regular expression
    //  Template partial specialization of findMatchingStrings
    template<class StringType>
    labelList findStrings
    (
        const std::string& rePattern,
        const UList<StringType>& lst,
        const bool invert=false
    )
    {
        const regExp re(rePattern);
        return findMatchingStrings(re, lst, invert);
    }

    //- Return list indices for strings matching the regular expression
    //  Template partial specialization of findMatchingStrings
    template<class StringType>
    labelList findStrings
    (
        const wordRe& wre,
        const UList<StringType>& lst,
        const bool invert=false
    )
    {
        return findMatchingStrings(wre, lst, invert);
    }


    //- Return list indices for strings matching one of the regular expression
    //  Template partial specialization of findMatchingStrings
    template<class StringType>
    labelList findStrings
    (
        const wordReListMatcher& matcher,
        const UList<StringType>& lst,
        const bool invert=false
    )
    {
        return findMatchingStrings(matcher, lst, invert);
    }

    // subsetting multi-string matches (similar to ListOp):

    //- Extract elements of StringList when regular expression matches
    //  optionally invert the match
    //  eg, to extract all selected elements:
    //    subsetMatchingStrings<regExp, stringList>(myRegExp, lst);
    template<class Matcher, class StringListType>
    StringListType subsetMatchingStrings
    (
        const Matcher&,
        const StringListType&,
        const bool invert=false
    );

    //- Extract elements of StringList when regular expression matches
    //  Template partial specialization of subsetMatchingStrings
    template<class StringListType>
    StringListType subsetStrings
    (
        const regExp& re,
        const StringListType& lst,
        const bool invert=false
    )
    {
        return subsetMatchingStrings(re, lst, invert);
    }

    //- Extract elements of StringList when regular expression matches
    //  Template partial specialization of subsetMatchingStrings
    template<class StringListType>
    StringListType subsetStrings
    (
        const char* rePattern,
        const StringListType& lst,
        const bool invert=false
    )
    {
        const regExp re(rePattern);
        return subsetMatchingStrings(re, lst, invert);
    }

    //- Extract elements of StringList when regular expression matches
    //  Template partial specialization of subsetMatchingStrings
    template<class StringListType>
    StringListType subsetStrings
    (
        const std::string& rePattern,
        const StringListType& lst,
        const bool invert=false
    )
    {
        const regExp re(rePattern);
        return subsetMatchingStrings(re, lst, invert);
    }

    //- Extract elements of StringList when regular expression matches
    //  Template partial specialization of subsetMatchingStrings
    template<class StringListType>
    StringListType subsetStrings
    (
        const wordRe& wre,
        const StringListType& lst,
        const bool invert=false
    )
    {
        return subsetMatchingStrings(wre, lst, invert);
    }

    //- Extract elements of StringList when regular expression matches
    //  Template partial specialization of subsetMatchingStrings
    template<class StringListType>
    StringListType subsetStrings
    (
        const wordReListMatcher& matcher,
        const StringListType& lst,
        const bool invert=false
    )
    {
        return subsetMatchingStrings(matcher, lst, invert);
    }


    //- Inplace extract elements of StringList when regular expression matches
    //  optionally invert the match
    //  eg, to extract all selected elements:
    //    inplaceSubsetMatchingStrings<regExp, stringList>(myRegExp, lst);
    template<class Matcher, class StringListType>
    void inplaceSubsetMatchingStrings
    (
        const Matcher&,
        StringListType&,
        const bool invert=false
    );

    //- Inplace extract elements of StringList when regular expression matches
    //  Template partial specialization of inplaceSubsetMatchingStrings
    template<class StringListType>
    void inplaceSubsetStrings
    (
        const regExp& re,
        StringListType& lst,
        const bool invert=false
    )
    {
        inplaceSubsetMatchingStrings(re, lst, invert);
    }

    //- Inplace extract elements of StringList when regular expression matches
    //  Template partial specialization of inplaceSubsetMatchingStrings
    template<class StringListType>
    void inplaceSubsetStrings
    (
        const char* rePattern,
        StringListType& lst,
        const bool invert=false
    )
    {
        const regExp re(rePattern);
        inplaceSubsetMatchingStrings(re, lst, invert);
    }

    //- Inplace extract elements of StringList when regular expression matches
    //  Template partial specialization of inplaceSubsetMatchingStrings
    template<class StringListType>
    void inplaceSubsetStrings
    (
        const std::string& rePattern,
        StringListType& lst,
        const bool invert=false
    )
    {
        const regExp re(rePattern);
        inplaceSubsetMatchingStrings(re, lst, invert);
    }

    //- Inplace extract elements of StringList when regular expression matches
    //  Template partial specialization of inplaceSubsetMatchingStrings
    template<class StringListType>
    void inplaceSubsetStrings
    (
        const wordRe& wre,
        StringListType& lst,
        const bool invert=false
    )
    {
        inplaceSubsetMatchingStrings(wre, lst, invert);
    }

    //- Inplace extract elements of StringList when regular expression matches
    //  Template partial specialization of inplaceSubsetMatchingStrings
    template<class StringListType>
    void inplaceSubsetStrings
    (
        const wordReListMatcher& matcher,
        StringListType& lst,
        const bool invert=false
    )
    {
        inplaceSubsetMatchingStrings(matcher, lst, invert);
    }

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "stringListOpsTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
