/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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

#include "LUscalarMatrix.H"
#include "SubField.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::LUscalarMatrix::solve
(
    Field<Type>& x,
    const Field<Type>& source
) const
{
    // If x and source are different initialise x = source
    if (&x != &source)
    {
        x = source;
    }

    if (Pstream::parRun())
    {
        Field<Type> X(m());

        if (Pstream::master(comm_))
        {
            typename Field<Type>::subField
            (
                X,
                x.size()
            ) = x;

            for
            (
                int slave=Pstream::firstSlave();
                slave<=Pstream::lastSlave(comm_);
                slave++
            )
            {
                IPstream::read
                (
                    Pstream::commsTypes::scheduled,
                    slave,
                    reinterpret_cast<char*>
                    (
                        &(X[procOffsets_[slave]])
                    ),
                    (procOffsets_[slave+1]-procOffsets_[slave])*sizeof(Type),
                    Pstream::msgType(),
                    comm_
                );
            }
        }
        else
        {
            OPstream::write
            (
                Pstream::commsTypes::scheduled,
                Pstream::masterNo(),
                reinterpret_cast<const char*>(x.begin()),
                x.byteSize(),
                Pstream::msgType(),
                comm_
            );
        }

        if (Pstream::master(comm_))
        {
            LUBacksubstitute(*this, pivotIndices_, X);

            x = typename Field<Type>::subField
            (
                X,
                x.size()
            );

            for
            (
                int slave=Pstream::firstSlave();
                slave<=Pstream::lastSlave(comm_);
                slave++
            )
            {
                OPstream::write
                (
                    Pstream::commsTypes::scheduled,
                    slave,
                    reinterpret_cast<const char*>
                    (
                        &(X[procOffsets_[slave]])
                    ),
                    (procOffsets_[slave + 1]-procOffsets_[slave])*sizeof(Type),
                    Pstream::msgType(),
                    comm_
                );
            }
        }
        else
        {
            IPstream::read
            (
                Pstream::commsTypes::scheduled,
                Pstream::masterNo(),
                reinterpret_cast<char*>(x.begin()),
                x.byteSize(),
                Pstream::msgType(),
                comm_
            );
        }
    }
    else
    {
        LUBacksubstitute(*this, pivotIndices_, x);
    }
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::LUscalarMatrix::solve
(
    const Field<Type>& source
) const
{
    tmp<Field<Type>> tx(new Field<Type>(m()));
    Field<Type>& x = tx.ref();

    solve(x, source);

    return tx;
}


// ************************************************************************* //
