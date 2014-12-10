#include "primitiveFields.H"
#include "Random.H"
#include "cpuTime.H"
#include "IOstreams.H"
#include "OFstream.H"

using namespace Foam;

int main()
{
    Info<< "Initialising fields" << endl;

    const label nIter = 100;
    const label size = 10000000;
    const label rndAddrSkip = 40;
    const label redFac = 6;
    const label redSize = size/redFac;
    Random genAddr(100);

    double* f1 = new double[size];
    double* f2 = new double[size];
    double* f3 = new double[size];
    double* f4 = new double[size];
    double* fr = new double[redSize];
    label*  addr = new label[size];
    label*  redAddr = new label[size];
    label*  redAddr2 = new label[size];

    for (register label i=0; i<size; i++)
    {
        f1[i] = 1.0;
        f2[i] = 1.0;
        f3[i] = 1.0;
        addr[i] = i;
        redAddr[i] = i/redFac;
        redAddr2[i] = (size - i - 1)/redFac;
    }

    for (register label i=0; i<size; i+=rndAddrSkip)
    {
        addr[i] = genAddr.integer(0, size-1);
    }

    for (register label i=0; i<redSize; i++)
    {
        fr[i] = 1.0;
    }

    Info<< "Done\n" <<endl;

    {
        Info<< "Single loop combined operation (expression templates)"
            << endl;

        cpuTime executionTime;

        for (int j=0; j<nIter; j++)
        {
            for (register label i=0; i<size; i++)
            {
                f4[i] = f1[i] + f2[i] - f3[i];
            }
        }

        Info<< "ExecutionTime = "
            << executionTime.elapsedCpuTime()
            << " s\n" << endl;

        Snull<< f4[1] << endl << endl;
    }

    {
        Info<< "Single loop combined operation with indirect addressing"
            << endl;

        cpuTime executionTime;

        for (int j=0; j<nIter; j++)
        {
            for (register label i=0; i<size; i++)
            {
                f4[addr[i]] = f1[addr[i]] + f2[addr[i]] - f3[addr[i]];
            }
        }

        Info<< "ExecutionTime = "
            << executionTime.elapsedCpuTime()
            << " s\n" << endl;

        Snull<< f4[1] << endl << endl;
    }

    {
        Info<< "Single loop reduction operation"
            << endl;

        cpuTime executionTime;
        label redOffset = (size - 1)/redFac;

        for (int j=0; j<nIter; j++)
        {
            for (register label i=0; i<size; i++)
            {
                label j = i/redFac;
                fr[j] += f1[i];
                fr[redOffset - j] -= f2[i];
            }
        }

        Info<< "ExecutionTime = "
            << executionTime.elapsedCpuTime()
            << " s\n" << endl;

        Snull<< fr[1] << endl << endl;
    }

    {
        Info<< "Single loop reduction operation with indirect addressing"
            << endl;

        cpuTime executionTime;

        for (int j=0; j<nIter; j++)
        {
            for (register label i=0; i<size; i++)
            {
                fr[redAddr[i]] += f1[i];
                fr[redAddr2[i]] -= f2[i];
            }
        }

        Info<< "ExecutionTime = "
            << executionTime.elapsedCpuTime()
            << " s\n" << endl;

        Snull<< fr[1] << endl << endl;
    }

    {
        Info<< "Separate loops ?= operations" << endl;

        cpuTime executionTime;

        for (int j=0; j<nIter; j++)
        {
            for (register label i=0; i<size; i++)
            {
                f4[i] = f1[i];
            }
            for (register label i=0; i<size; i++)
            {
                f4[i] += f2[i];
            }
            for (register label i=0; i<size; i++)
            {
                f4[i] -= f3[i];
            }
        }

        Info<< "ExecutionTime = "
            << executionTime.elapsedCpuTime()
            << " s\n" << endl;

        Snull<< f4[1] << endl << endl;
    }

    {
        Info<< "OpenFOAM field algebra" << endl;

        scalarField
            sf1(size, 1.0),
            sf2(size, 1.0),
            sf3(size, 1.0),
            sf4(size);

        cpuTime executionTime;

        for (int j=0; j<nIter; j++)
        {
            //sf4 = sf1 + sf2 - sf3;
            sf4 = sf1;
            sf4 += sf2;
            sf4 -= sf3;
        }

        Info<< "ExecutionTime = "
            << executionTime.elapsedCpuTime()
            << " s\n" << endl;

        Snull<< sf4[1] << endl << endl;
    }
}
