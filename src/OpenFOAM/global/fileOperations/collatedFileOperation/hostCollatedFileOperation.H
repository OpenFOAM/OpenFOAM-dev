/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017-2018 OpenFOAM Foundation
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
    Foam::fileOperations::hostCollatedFileOperation

Description
    Version of collatedFileOperation with multiple read/write ranks.

    In parallel it will assume ranks are sorted according to hostname
    and the lowest rank per hostname will be the IO rank. The output directories
    will get a unique name processors<N>_<low>-<high> where N is the overall
    number of processors and low and high is the range of ranks contained
    in the files. Each of these subsets uses its own communicator.

    Instead of using the hostnames the IO ranks can be assigned using the
    FOAM_IORANKS environment variable (also when running non-parallel), e.g.
    when decomposing into 4:

        FOAM_IORANKS='(0 2)' decomposePar -fileHandler hostCollated

    will generate

        processors4_0-1/
            containing data for processors 0 to 1
        processors4_2-3/
            containing data for processors 2 to 3

See also
    collatedFileOperation

SourceFiles
    hostCollatedFileOperation.C

\*---------------------------------------------------------------------------*/

#ifndef fileOperations_hostCollatedFileOperation_H
#define fileOperations_hostCollatedFileOperation_H

#include "collatedFileOperation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fileOperations
{

/*---------------------------------------------------------------------------*\
                 Class hostCollatedFileOperation Declaration
\*---------------------------------------------------------------------------*/

class hostCollatedFileOperation
:
    public collatedFileOperation
{
   // Private Member Functions

        //- Get the list of processors part of this set
        static labelList subRanks(const label n);


public:

        //- Runtime type information
        TypeName("hostCollated");


    // Constructors

        //- Construct null
        hostCollatedFileOperation(const bool verbose);


    //- Destructor
    virtual ~hostCollatedFileOperation();
};


/*---------------------------------------------------------------------------*\
            Class hostCollatedFileOperationInitialise Declaration
\*---------------------------------------------------------------------------*/

class hostCollatedFileOperationInitialise
:
    public collatedFileOperationInitialise
{
public:

    // Constructors

        //- Construct from components
        hostCollatedFileOperationInitialise(int& argc, char**& argv)
        :
            collatedFileOperationInitialise(argc, argv)
        {}


    //- Destructor
    virtual ~hostCollatedFileOperationInitialise()
    {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fileOperations
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
