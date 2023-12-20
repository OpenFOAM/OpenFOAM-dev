/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2023 OpenFOAM Foundation
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

Class
    Foam::nonConformalPolyPatch

Description
    Non-conformal poly patch. This patch is a placeholder and must have no
    faces. This patch is linked to an "original" patch. Other parts of the code
    (i.e., finite volume) will use this patch to create non-empty patches that
    do not conform to the poly mesh.

SourceFiles
    nonConformalPolyPatch.C

\*---------------------------------------------------------------------------*/

#ifndef nonConformalPolyPatch_H
#define nonConformalPolyPatch_H

#include "cyclicPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class nonConformalPolyPatch Declaration
\*---------------------------------------------------------------------------*/

class nonConformalPolyPatch
{
private:

    // Private

        //- Reference to the polyPatch
        const polyPatch& patch_;

        //- Name of the patch from which this patch originates
        mutable word origPatchName_;

        //- Index of the patch from which this patch originates
        mutable label origPatchIndex_;


    // Private Member Functions

        //- Validate the size of the patch
        void validateSize() const;


protected:

    // Protected Member Functions

        //- Reset the patch name
        virtual void rename(const wordList& newNames);

        //- Reset the patch index
        virtual void reorder(const labelUList& newToOldIndex);


public:

    //- Runtime type information
    TypeName("nonConformal");


    // Constructors

        //- Construct from a patch
        nonConformalPolyPatch(const polyPatch& patch);

        //- Construct from a patch and an original patch name
        nonConformalPolyPatch
        (
            const polyPatch& patch,
            const word& origPatchName
        );

        //- Construct from a patch and a dictionary
        nonConformalPolyPatch
        (
            const polyPatch& patch,
            const dictionary& dict
        );

        //- Construct from a patch and a non-conformal patch
        nonConformalPolyPatch
        (
            const polyPatch& patch,
            const nonConformalPolyPatch& ncPatch
        );


    //- Destructor
    virtual ~nonConformalPolyPatch();


    // Member Functions

        //- Reference to the polyPatch
        const polyPatch& patch() const;

        //- Original patch name
        const word& origPatchName() const;

        //- Original patchID
        label origPatchIndex() const;

        //- Original patch
        const polyPatch& origPatch() const;

        //- Write the polyPatch data as a dictionary
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
