/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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
    Foam::multiRegionPrefixer

SourceFiles
    multiRegionPrefixer.C

\*---------------------------------------------------------------------------*/

#ifndef multiRegionPrefixer_H
#define multiRegionPrefixer_H

#include "wordList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class multiRegionPrefixer Declaration
\*---------------------------------------------------------------------------*/

class multiRegionPrefixer
{
private:

    // Private Data

        //- Should we prefix a single region?
        const bool prefixSingleRegion_;

        //- Region names
        const wordList& regionNames_;


    // Private Member Functions

        //- Are we prefixing?
        bool prefixes() const;

        //- Width of the write prefix
        string::size_type prefixWidth() const;


public:

    // Public classes

        //- Prefixer for a single region
        class regionPrefixer
        {
        private:

            // Private data

                //- The multi-region prefixing engine
                const multiRegionPrefixer& mrp_;

                //- The region index
                label regioni_;


        public:

            // Constructors

                //- Construct from components
                regionPrefixer
                (
                    const multiRegionPrefixer& mrp,
                    const label regioni
                );

                //- Disallow copy construct
                regionPrefixer(const regionPrefixer& rp) = delete;

                //- Move construct
                regionPrefixer(regionPrefixer&& rp);


            //- Destructor
            ~regionPrefixer();


            // Member operators

                //- Disallow copy assign
                regionPrefixer& operator=(const regionPrefixer&) = delete;

                //- Disallow move assign
                regionPrefixer& operator=(regionPrefixer&&) = delete;
        };


    // Constructors

        //- Construct from components
        multiRegionPrefixer
        (
            const bool prefixSingleRegion,
            const wordList& regionNames
        );


    //- Destructor
    ~multiRegionPrefixer();
};


/*---------------------------------------------------------------------------*\
                         Class RegionConstRef Declaration
\*---------------------------------------------------------------------------*/

template<class Region>
class RegionConstRef
:
    private multiRegionPrefixer::regionPrefixer
{
private:

    // Private Data

        //- Reference to the region object
        const Region& r_;


public:

    // Constructors

        //- Construct form components
        inline RegionConstRef
        (
            const multiRegionPrefixer& mrp,
            const label regioni,
            const Region& r
        )
        :
            multiRegionPrefixer::regionPrefixer(mrp, regioni),
            r_(r)
        {}


    // Member operators

        //- Cast to reference
        inline operator const Region&() const
        {
            return r_;
        }

        //- Obtain the reference
        inline const Region& operator()() const
        {
            return r_;
        }
};


/*---------------------------------------------------------------------------*\
                         Class RegionRef Declaration
\*---------------------------------------------------------------------------*/

template<class Region>
class RegionRef
:
    public RegionConstRef<Region>
{
private:

    // Private Data

        //- Reference to the region object
        Region& r_;


public:

    // Constructors

        //- Construct form components
        inline RegionRef
        (
            const multiRegionPrefixer& mrp,
            const label regioni,
            Region& r
        )
        :
            RegionConstRef<Region>(mrp, regioni, r),
            r_(r)
        {}


    // Member operators

        //- Cast to reference
        inline operator Region&() const
        {
            return r_;
        }

        //- Obtain the reference
        inline Region& operator()()
        {
            return r_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
