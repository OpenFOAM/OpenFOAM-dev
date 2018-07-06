/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "argList.H"
#include "OSspecific.H"
#include "clock.H"
#include "IFstream.H"
#include "dictionary.H"
#include "IOobject.H"
#include "jobInfo.H"
#include "labelList.H"
#include "regIOobject.H"
#include "dynamicCode.H"
#include "fileOperation.H"
#include "fileOperationInitialise.H"
#include "stringListOps.H"

#include <cctype>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

Foam::SLList<Foam::string>    Foam::argList::validArgs;
Foam::HashTable<Foam::string> Foam::argList::validOptions;
Foam::HashTable<Foam::string> Foam::argList::validParOptions;
Foam::HashTable<Foam::string> Foam::argList::optionUsage;
Foam::SLList<Foam::string>    Foam::argList::notes;
Foam::string::size_type Foam::argList::usageMin = 20;
Foam::string::size_type Foam::argList::usageMax = 80;
Foam::word Foam::argList::postProcessOptionName("postProcess");

Foam::argList::initValidTables::initValidTables()
{
    argList::addOption
    (
        "case", "dir",
        "specify alternate case directory, default is the cwd"
    );
    argList::addBoolOption("parallel", "run in parallel");
    validParOptions.set("parallel", "");
    argList::addOption
    (
        "roots", "(dir1 .. dirN)",
        "slave root directories for distributed running"
    );
    validParOptions.set("roots", "(dir1 .. dirN)");

    argList::addOption
    (
        "hostRoots", "(((host1 dir1) .. (hostN dirN))",
        "slave root directories (per host) for distributed running"
    );
    validParOptions.set("hostRoots", "((host1 dir1) .. (hostN dirN))");

    argList::addBoolOption
    (
        "noFunctionObjects",
        "do not execute functionObjects"
    );

    argList::addOption
    (
        "fileHandler",
        "handler",
        "override the fileHandler"
     );

    Pstream::addValidParOptions(validParOptions);
}


void Foam::argList::initValidTables::clear()
{
    argList::removeOption("case");
    argList::removeOption("parallel");
    argList::removeOption("roots");
    argList::removeOption("hostRoots");
    argList::removeOption("noFunctionObjects");
    argList::removeOption("fileHandler");
}


Foam::argList::initValidTables dummyInitValidTables;


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

void Foam::argList::addBoolOption
(
    const word& opt,
    const string& usage
)
{
    addOption(opt, "", usage);
}


void Foam::argList::addOption
(
    const word& opt,
    const string& param,
    const string& usage
)
{
    validOptions.set(opt, param);
    if (!usage.empty())
    {
        optionUsage.set(opt, usage);
    }
}


void Foam::argList::addUsage
(
    const word& opt,
    const string& usage
)
{
    if (usage.empty())
    {
        optionUsage.erase(opt);
    }
    else
    {
        optionUsage.set(opt, usage);
    }
}


void Foam::argList::addNote(const string& note)
{
    if (!note.empty())
    {
        notes.append(note);
    }
}


void Foam::argList::removeOption(const word& opt)
{
    validOptions.erase(opt);
    optionUsage.erase(opt);
}


void Foam::argList::noParallel()
{
    removeOption("parallel");
    removeOption("roots");
    removeOption("hostRoots");
    validParOptions.clear();
}


void Foam::argList::printOptionUsage
(
    const label location,
    const string& str
)
{
    const string::size_type textWidth = usageMax - usageMin;
    const string::size_type strLen = str.size();

    if (strLen)
    {
        // Minimum of 2 spaces between option and usage:
        if (string::size_type(location) + 2 <= usageMin)
        {
            for (string::size_type i = location; i < usageMin; ++i)
            {
                Info<<' ';
            }
        }
        else
        {
            // or start a new line
            Info<< nl;
            for (string::size_type i = 0; i < usageMin; ++i)
            {
                Info<<' ';
            }
        }

        // Text wrap
        string::size_type pos = 0;
        while (pos != string::npos && pos + textWidth < strLen)
        {
            // Potential end point and next point
            string::size_type curr = pos + textWidth - 1;
            string::size_type next = string::npos;

            if (isspace(str[curr]))
            {
                // We were lucky: ended on a space
                next = str.find_first_not_of(" \t\n", curr);
            }
            else if (isspace(str[curr+1]))
            {
                // The next one is a space - so we are okay
                curr++;  // otherwise the length is wrong
                next = str.find_first_not_of(" \t\n", curr);
            }
            else
            {
                // Search for end of a previous word break
                string::size_type prev = str.find_last_of(" \t\n", curr);

                // Reposition to the end of previous word if possible
                if (prev != string::npos && prev > pos)
                {
                    curr = prev;
                }
            }

            if (next == string::npos)
            {
                next = curr + 1;
            }

            // Indent following lines (not the first one)
            if (pos)
            {
                for (string::size_type i = 0; i < usageMin; ++i)
                {
                    Info<<' ';
                }
            }

            Info<< str.substr(pos, (curr - pos)).c_str() << nl;
            pos = next;
        }

        // Output the remainder of the string
        if (pos != string::npos)
        {
            // Indent following lines (not the first one)
            if (pos)
            {
                for (string::size_type i = 0; i < usageMin; ++i)
                {
                    Info<<' ';
                }
            }

            Info<< str.substr(pos).c_str() << nl;
        }
    }
    else
    {
        Info<< nl;
    }
}


bool Foam::argList::postProcess(int argc, char *argv[])
{
    bool postProcessOption = false;

    for (int i=1; i<argc; i++)
    {
        postProcessOption = argv[i] == '-' + postProcessOptionName;
        if (postProcessOption) break;
    }

    return postProcessOption;
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Convert argv -> args_
// Transform sequences with "(" ... ")" into string lists in the process
bool Foam::argList::regroupArgv(int& argc, char**& argv)
{
    int nArgs = 0;
    int listDepth = 0;
    string tmpString;

    // Note: we also re-write directly into args_
    // and use a second pass to sort out args/options
    for (int argI = 0; argI < argc; ++argI)
    {
        if (strcmp(argv[argI], "(") == 0)
        {
            ++listDepth;
            tmpString += "(";
        }
        else if (strcmp(argv[argI], ")") == 0)
        {
            if (listDepth)
            {
                --listDepth;
                tmpString += ")";
                if (listDepth == 0)
                {
                    args_[nArgs++] = tmpString;
                    tmpString.clear();
                }
            }
            else
            {
                args_[nArgs++] = argv[argI];
            }
        }
        else if (listDepth)
        {
            // Quote each string element
            tmpString += "\"";
            tmpString += argv[argI];
            tmpString += "\"";
        }
        else
        {
            args_[nArgs++] = argv[argI];
        }
    }

    if (tmpString.size())
    {
        args_[nArgs++] = tmpString;
    }

    args_.setSize(nArgs);

    return nArgs < argc;
}


void Foam::argList::getRootCase()
{
    fileName casePath;

    // [-case dir] specified
    HashTable<string>::const_iterator iter = options_.find("case");

    if (iter != options_.end())
    {
        casePath = iter();
        casePath.clean();

        if (casePath.empty() || casePath == ".")
        {
            // Handle degenerate form and '-case .' like no -case specified
            casePath = cwd();
            options_.erase("case");
        }
        else if (!casePath.isAbsolute() && casePath.name() == "..")
        {
            // Avoid relative cases ending in '..' - makes for very ugly names
            casePath = cwd()/casePath;
            casePath.clean();
        }
    }
    else
    {
        // Nothing specified, use the current dir
        casePath = cwd();
    }

    rootPath_   = casePath.path();
    globalCase_ = casePath.name();
    case_       = globalCase_;


    // Set the case and case-name as an environment variable
    if (rootPath_.isAbsolute())
    {
        // Absolute path - use as-is
        setEnv("FOAM_CASE", rootPath_/globalCase_, true);
        setEnv("FOAM_CASENAME", globalCase_, true);
    }
    else
    {
        // Qualify relative path
        casePath = cwd()/rootPath_/globalCase_;
        casePath.clean();

        setEnv("FOAM_CASE", casePath, true);
        setEnv("FOAM_CASENAME", casePath.name(), true);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::argList::argList
(
    int& argc,
    char**& argv,
    bool checkArgs,
    bool checkOpts,
    const bool initialise
)
:
    args_(argc),
    options_(argc)
{
    // Check for fileHandler
    word handlerType(getEnv("FOAM_FILEHANDLER"));
    for (int argI = 0; argI < argc; ++argI)
    {
        if (argv[argI][0] == '-')
        {
            const char *optionName = &argv[argI][1];
            if (string(optionName) == "fileHandler")
            {
                handlerType = argv[argI+1];
                break;
            }
        }
    }
    if (handlerType.empty())
    {
        handlerType = fileOperation::defaultFileHandler;
    }

    // Detect any parallel options
    bool needsThread = fileOperations::fileOperationInitialise::New
    (
        handlerType,
        argc,
        argv
    )().needsThreading();


    // Check if this run is a parallel run by searching for any parallel option
    // If found call runPar which might filter argv
    for (int argI = 0; argI < argc; ++argI)
    {
        if (argv[argI][0] == '-')
        {
            const char *optionName = &argv[argI][1];

            if (validParOptions.found(optionName))
            {
                parRunControl_.runPar(argc, argv, needsThread);
                break;
            }
        }
    }

    // Convert argv -> args_ and capture ( ... ) lists
    // for normal arguments and for options
    regroupArgv(argc, argv);

    // Get executable name
    args_[0]    = fileName(argv[0]);
    executable_ = fileName(argv[0]).name();

    // Check arguments and options, we already have argv[0]
    int nArgs = 1;
    argListStr_ = args_[0];

    for (int argI = 1; argI < args_.size(); ++argI)
    {
        argListStr_ += ' ';
        argListStr_ += args_[argI];

        if (args_[argI][0] == '-')
        {
            const char *optionName = &args_[argI][1];

            if
            (
                (
                    validOptions.found(optionName)
                 && !validOptions[optionName].empty()
                )
             || (
                    validParOptions.found(optionName)
                 && !validParOptions[optionName].empty()
                )
            )
            {
                ++argI;
                if (argI >= args_.size())
                {
                    FatalError
                        <<"Option '-" << optionName
                        << "' requires an argument" << endl;
                    printUsage();
                    FatalError.exit();
                }

                argListStr_ += ' ';
                argListStr_ += args_[argI];
                options_.insert(optionName, args_[argI]);
            }
            else
            {
                options_.insert(optionName, "");
            }
        }
        else
        {
            if (nArgs != argI)
            {
                args_[nArgs] = args_[argI];
            }
            ++nArgs;
        }
    }

    args_.setSize(nArgs);

    parse(checkArgs, checkOpts, initialise);
}


Foam::argList::argList
(
    const argList& args,
    const HashTable<string>& options,
    bool checkArgs,
    bool checkOpts,
    bool initialise
)
:
    parRunControl_(args.parRunControl_),
    args_(args.args_),
    options_(options),
    executable_(args.executable_),
    rootPath_(args.rootPath_),
    globalCase_(args.globalCase_),
    case_(args.case_),
    argListStr_(args.argListStr_)
{
    parse(checkArgs, checkOpts, initialise);
}


void Foam::argList::parse
(
    bool checkArgs,
    bool checkOpts,
    bool initialise
)
{
    // Help/documentation options:
    //   -help    print the usage
    //   -doc     display application documentation in browser
    //   -srcDoc  display source code in browser
    if
    (
        options_.found("help")
     || options_.found("doc")
     || options_.found("srcDoc")
    )
    {
        if (options_.found("help"))
        {
            printUsage();
        }

        // Only display one or the other
        if (options_.found("srcDoc"))
        {
            displayDoc(true);
        }
        else if (options_.found("doc"))
        {
            displayDoc(false);
        }

        ::exit(0);
    }

    // Print the usage message and exit if the number of arguments is incorrect
    if (!check(checkArgs, checkOpts))
    {
        FatalError.exit();
    }


    if (initialise)
    {
        string dateString = clock::date();
        string timeString = clock::clockTime();

        // Print the banner once only for parallel runs
        if (Pstream::master() && writeInfoHeader)
        {
            IOobject::writeBanner(Info, true)
                << "Build  : " << Foam::FOAMbuild << nl
                << "Exec   : " << argListStr_.c_str() << nl
                << "Date   : " << dateString.c_str() << nl
                << "Time   : " << timeString.c_str() << nl
                << "Host   : " << hostName() << nl
                << "PID    : " << pid() << endl;
        }

        jobInfo_.add("startDate", dateString);
        jobInfo_.add("startTime", timeString);
        jobInfo_.add("userName", userName());
        jobInfo_.add("foamVersion", word(FOAMversion));
        jobInfo_.add("code", executable_);
        jobInfo_.add("argList", argListStr_);
        jobInfo_.add("currentDir", cwd());
        jobInfo_.add("PPID", ppid());
        jobInfo_.add("PGID", pgid());

        // Add build information - only use the first word
        {
            std::string build(Foam::FOAMbuild);
            std::string::size_type found = build.find(' ');
            if (found != std::string::npos)
            {
                build.resize(found);
            }
            jobInfo_.add("foamBuild", build);
        }
    }


    // Set fileHandler. In increasing order of priority:
    // 1. default = uncollated
    // 2. environment var FOAM_FILEHANDLER
    // 3. etc/controlDict optimisationSwitches 'fileHandler'
    // 4. system/controlDict 'fileHandler' (not handled here; done in TimeIO.C)

    {
        word handlerType(getEnv("FOAM_FILEHANDLER"));
        HashTable<string>::const_iterator iter = options_.find("fileHandler");
        if (iter != options_.end())
        {
            handlerType = iter();
        }

        if (handlerType.empty())
        {
            handlerType = fileOperation::defaultFileHandler;
        }

        autoPtr<fileOperation> handler
        (
            fileOperation::New
            (
                handlerType,
                writeInfoHeader
            )
        );
        Foam::fileHandler(handler);
    }


    stringList slaveMachine;
    stringList slaveProcs;

    // Collect slave machine/pid
    if (parRunControl_.parRun())
    {
        if (Pstream::master())
        {
            slaveMachine.setSize(Pstream::nProcs() - 1);
            slaveProcs.setSize(Pstream::nProcs() - 1);
            label proci = 0;
            for
            (
                int slave = Pstream::firstSlave();
                slave <= Pstream::lastSlave();
                slave++
            )
            {
                IPstream fromSlave(Pstream::commsTypes::scheduled, slave);

                string slaveBuild;
                label slavePid;
                fromSlave >> slaveBuild >> slaveMachine[proci] >> slavePid;

                slaveProcs[proci] = slaveMachine[proci]+"."+name(slavePid);
                proci++;

                // Check build string to make sure all processors are running
                // the same build
                if (slaveBuild != Foam::FOAMbuild)
                {
                    FatalErrorIn(executable())
                        << "Master is running version " << Foam::FOAMbuild
                        << "; slave " << proci << " is running version "
                        << slaveBuild
                        << exit(FatalError);
                }
            }
        }
        else
        {
            OPstream toMaster
            (
                Pstream::commsTypes::scheduled,
                Pstream::masterNo()
            );
            toMaster << string(Foam::FOAMbuild) << hostName() << pid();
        }
    }


    // Case is a single processor run unless it is running parallel
    int nProcs = 1;

    // Roots if running distributed
    fileNameList roots;

    // If this actually is a parallel run
    if (parRunControl_.parRun())
    {
        // For the master
        if (Pstream::master())
        {
            // Establish rootPath_/globalCase_/case_ for master
            getRootCase();

            // See if running distributed (different roots for different procs)
            label dictNProcs = -1;
            fileName source;

            if (options_.found("roots"))
            {
                source = "-roots";
                IStringStream is(options_["roots"]);
                roots = readList<fileName>(is);

                if (roots.size() != 1)
                {
                    dictNProcs = roots.size()+1;
                }
            }
            else if (options_.found("hostRoots"))
            {
                source = "-hostRoots";
                IStringStream is(options_["hostRoots"]);
                List<Tuple2<wordRe, fileName>> hostRoots(is);

                roots.setSize(Pstream::nProcs()-1);
                forAll(hostRoots, i)
                {
                    const Tuple2<wordRe, fileName>& hostRoot = hostRoots[i];
                    const wordRe& re = hostRoot.first();
                    labelList matchedRoots(findStrings(re, slaveMachine));
                    forAll(matchedRoots, matchi)
                    {
                        label slavei = matchedRoots[matchi];
                        if (roots[slavei] != wordRe())
                        {
                            FatalErrorInFunction
                                << "Slave " << slaveMachine[slavei]
                                << " has multiple matching roots in "
                                << hostRoots << exit(FatalError);
                        }
                        else
                        {
                            roots[slavei] = hostRoot.second();
                        }
                    }
                }

                // Check
                forAll(roots, slavei)
                {
                    if (roots[slavei] == wordRe())
                    {
                        FatalErrorInFunction
                            << "Slave " << slaveMachine[slavei]
                            << " has no matching roots in "
                            << hostRoots << exit(FatalError);
                    }
                }

                if (roots.size() != 1)
                {
                    dictNProcs = roots.size()+1;
                }
            }
            else
            {
                source = rootPath_/globalCase_/"system/decomposeParDict";
                IFstream decompDictStream(source);

                if (!decompDictStream.good())
                {
                    FatalError
                        << "Cannot read "
                        << decompDictStream.name()
                        << exit(FatalError);
                }

                dictionary decompDict(decompDictStream);

                dictNProcs = readLabel
                (
                    decompDict.lookup("numberOfSubdomains")
                );

                if (decompDict.lookupOrDefault("distributed", false))
                {
                    decompDict.lookup("roots") >> roots;
                }
            }

            // Convenience:
            // when a single root is specified, use it for all processes
            if (roots.size() == 1)
            {
                const fileName rootName(roots[0]);
                roots.setSize(Pstream::nProcs()-1, rootName);

                // adjust dictNProcs for command-line '-roots' option
                if (dictNProcs < 0)
                {
                    dictNProcs = roots.size()+1;
                }
            }


            // Check number of processors.
            // nProcs     => number of actual procs
            // dictNProcs => number of procs specified in decompositionDict
            // nProcDirs  => number of processor directories
            //               (n/a when running distributed)
            //
            // - normal running : nProcs = dictNProcs = nProcDirs
            // - decomposition to more  processors : nProcs = dictNProcs
            // - decomposition to fewer processors : nProcs = nProcDirs
            if (dictNProcs > Pstream::nProcs())
            {
                FatalError
                    << source
                    << " specifies " << dictNProcs
                    << " processors but job was started with "
                    << Pstream::nProcs() << " processors."
                    << exit(FatalError);
            }


            // Distributed data
            if (roots.size())
            {
                if (roots.size() != Pstream::nProcs()-1)
                {
                    FatalError
                        << "number of entries in roots "
                        << roots.size()
                        << " is not equal to the number of slaves "
                        << Pstream::nProcs()-1
                        << exit(FatalError);
                }

                forAll(roots, i)
                {
                    roots[i].expand();
                }

                // Distribute the master's argument list (with new root)
                bool hadCaseOpt = options_.found("case");
                for
                (
                    int slave = Pstream::firstSlave();
                    slave <= Pstream::lastSlave();
                    slave++
                )
                {
                    options_.set("case", roots[slave-1]/globalCase_);

                    OPstream toSlave(Pstream::commsTypes::scheduled, slave);
                    toSlave << args_ << options_;
                }
                options_.erase("case");

                // Restore [-case dir]
                if (hadCaseOpt)
                {
                    options_.set("case", rootPath_/globalCase_);
                }
            }
            else
            {
                // Possibly going to fewer processors.
                // Check if all procDirs are there.
                if (dictNProcs < Pstream::nProcs())
                {
                    label nProcDirs = 0;
                    while
                    (
                        isDir
                        (
                            rootPath_/globalCase_/"processor"
                          + name(++nProcDirs)
                        )
                    )
                    {}

                    if (nProcDirs != Pstream::nProcs())
                    {
                        FatalError
                            << "number of processor directories = "
                            << nProcDirs
                            << " is not equal to the number of processors = "
                            << Pstream::nProcs()
                            << exit(FatalError);
                    }
                }

                // Distribute the master's argument list (unaltered)
                for
                (
                    int slave = Pstream::firstSlave();
                    slave <= Pstream::lastSlave();
                    slave++
                )
                {
                    OPstream toSlave(Pstream::commsTypes::scheduled, slave);
                    toSlave << args_ << options_;
                }
            }
        }
        else
        {
            // Collect the master's argument list
            IPstream fromMaster
            (
                Pstream::commsTypes::scheduled,
                Pstream::masterNo()
            );
            fromMaster >> args_ >> options_;

            // Establish rootPath_/globalCase_/case_ for slave
            getRootCase();
        }

        nProcs = Pstream::nProcs();
        case_ = globalCase_/(word("processor") + name(Pstream::myProcNo()));
    }
    else
    {
        // Establish rootPath_/globalCase_/case_
        getRootCase();
        case_ = globalCase_;
    }

    if (Pstream::master() && writeInfoHeader)
    {
        Info<< "Case   : " << (rootPath_/globalCase_).c_str() << nl
            << "nProcs : " << nProcs << endl;

        if (parRunControl_.parRun())
        {
            Info<< "Slaves : " << slaveProcs << nl;
            if (roots.size())
            {
                Info<< "Roots  : " << roots << nl;
            }
            Info<< "Pstream initialized with:" << nl
                << "    floatTransfer      : " << Pstream::floatTransfer << nl
                << "    nProcsSimpleSum    : " << Pstream::nProcsSimpleSum << nl
                << "    commsType          : "
                << Pstream::commsTypeNames[Pstream::defaultCommsType] << nl
                << "    polling iterations : " << Pstream::nPollProcInterfaces
                << endl;
        }
    }

    if (initialise)
    {
        jobInfo_.add("root", rootPath_);
        jobInfo_.add("case", globalCase_);
        jobInfo_.add("nProcs", nProcs);
        if (slaveProcs.size())
        {
            jobInfo_.add("slaves", slaveProcs);
        }
        if (roots.size())
        {
            jobInfo_.add("roots", roots);
        }
        jobInfo_.write(executable_, rootPath_/globalCase_);

        // Switch on signal trapping. We have to wait until after Pstream::init
        // since this sets up its own ones.
        sigFpe_.set(writeInfoHeader);
        sigInt_.set(writeInfoHeader);
        sigQuit_.set(writeInfoHeader);
        sigSegv_.set(writeInfoHeader);

        if (writeInfoHeader)
        {
            Info<< "fileModificationChecking : "
                << "Monitoring run-time modified files using "
                << regIOobject::fileCheckTypesNames
                    [
                        regIOobject::fileModificationChecking
                    ];
            if
            (
                (
                    regIOobject::fileModificationChecking
                 == regIOobject::timeStamp
                )
             || (
                    regIOobject::fileModificationChecking
                 == regIOobject::timeStampMaster
                )
            )
            {
                Info<< " (fileModificationSkew "
                    << regIOobject::fileModificationSkew << ")";
            }
            Info<< endl;

            Info<< "allowSystemOperations : ";
            if (dynamicCode::allowSystemOperations)
            {
                Info<< "Allowing user-supplied system call operations" << endl;
            }
            else
            {
                Info<< "Disallowing user-supplied system call operations"
                    << endl;
            }
        }

        if (Pstream::master() && writeInfoHeader)
        {
            Info<< endl;
            IOobject::writeDivider(Info);
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::argList::~argList()
{
    jobInfo_.end();

    // Delete file handler to flush any remaining IO
    autoPtr<fileOperation> dummy(nullptr);
    fileHandler(dummy);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::argList::setOption(const word& opt, const string& param)
{
    bool changed = false;

    // Only allow valid options
    if (validOptions.found(opt))
    {
        // Some options are to be protected
        if
        (
            opt == "case"
         || opt == "parallel"
         || opt == "roots"
         || opt == "hostRoots"
        )
        {
            FatalError
                <<"used argList::setOption on a protected option: '"
                << opt << "'" << endl;
            FatalError.exit();
        }

        if (validOptions[opt].empty())
        {
            // Bool option
            if (!param.empty())
            {
                // Disallow change of type
                FatalError
                    <<"used argList::setOption to change bool to non-bool: '"
                    << opt << "'" << endl;
                FatalError.exit();
            }
            else
            {
                // Did not previously exist
                changed = !options_.found(opt);
            }
        }
        else
        {
            // Non-bool option
            if (param.empty())
            {
                // Disallow change of type
                FatalError
                    <<"used argList::setOption to change non-bool to bool: '"
                    << opt << "'" << endl;
                FatalError.exit();
            }
            else
            {
                // Existing value needs changing, or did not previously exist
                changed = options_.found(opt) ? options_[opt] != param : true;
            }
        }
    }
    else
    {
        FatalError
            <<"used argList::setOption on an invalid option: '"
            << opt << "'" << nl << "allowed are the following:"
            << validOptions << endl;
        FatalError.exit();
    }

    // Set/change the option as required
    if (changed)
    {
        options_.set(opt, param);
    }

    return changed;
}


bool Foam::argList::unsetOption(const word& opt)
{
    // Only allow valid options
    if (validOptions.found(opt))
    {
        // Some options are to be protected
        if
        (
            opt == "case"
         || opt == "parallel"
         || opt == "roots"
         || opt == "hostRoots"
        )
        {
            FatalError
                <<"used argList::unsetOption on a protected option: '"
                << opt << "'" << endl;
            FatalError.exit();
        }

        // Remove the option, return true if state changed
        return options_.erase(opt);
    }
    else
    {
        FatalError
            <<"used argList::unsetOption on an invalid option: '"
            << opt << "'" << nl << "allowed are the following:"
            << validOptions << endl;
        FatalError.exit();
    }

    return false;
}


void Foam::argList::printNotes() const
{
    // Output notes directly - no automatic text wrapping
    if (!notes.empty())
    {
        Info<< nl;
        forAllConstIter(SLList<string>, notes, iter)
        {
            Info<< iter().c_str() << nl;
        }
    }
}


void Foam::argList::printUsage() const
{
    Info<< "\nUsage: " << executable_ << " [OPTIONS]";

    forAllConstIter(SLList<string>, validArgs, iter)
    {
        Info<< " <" << iter().c_str() << '>';
    }

    Info<< "\noptions:\n";

    wordList opts = validOptions.sortedToc();
    forAll(opts, optI)
    {
        const word& optionName = opts[optI];

        HashTable<string>::const_iterator iter = validOptions.find(optionName);
        Info<< "  -" << optionName;
        label len = optionName.size() + 3;  // Length includes leading '  -'

        if (iter().size())
        {
            // Length includes space and between option/param and '<>'
            len += iter().size() + 3;
            Info<< " <" << iter().c_str() << '>';
        }

        HashTable<string>::const_iterator usageIter =
            optionUsage.find(optionName);

        if (usageIter != optionUsage.end())
        {
            printOptionUsage
            (
                len,
                usageIter()
            );
        }
        else
        {
            Info<< nl;
        }
    }

    // Place srcDoc/doc/help options at the end
    Info<< "  -srcDoc";
    printOptionUsage
    (
        9,
        "display source code in browser"
    );

    Info<< "  -doc";
    printOptionUsage
    (
        6,
        "display application documentation in browser"
    );

    Info<< "  -help";
    printOptionUsage
    (
        7,
        "print the usage"
    );


    printNotes();

    Info<< nl
        <<"Using: OpenFOAM-" << Foam::FOAMversion
        << " (see www.OpenFOAM.org)" << nl
        <<"Build: " << Foam::FOAMbuild << nl
        << endl;
}


void Foam::argList::displayDoc(bool source) const
{
    const dictionary& docDict = debug::controlDict().subDict("Documentation");
    List<fileName> docDirs(docDict.lookup("doxyDocDirs"));
    fileName docExt(docDict.lookup("doxySourceFileExt"));

    // For source code: change foo_8C.html to foo_8C_source.html
    if (source)
    {
        docExt.replace(".", "_source.");
    }

    fileName docFile;
    fileName httpServer;
    bool found = false;

    forAll(docDirs, dirI)
    {
        // An HTTP server is treated as a special case ...
        if (docDirs[dirI].component(0) == "http:")
        {
            httpServer = docDirs[dirI]/executable_ + docExt;
        }
        else
        {
            // ... all other entries are treated as local directories

            // Remove the optional "file://"
            if (docDirs[dirI].component(0) == "file:")
            {
                docDirs[dirI].replace("file://", string::null);
            }


            // Expand the file name
            docFile = docDirs[dirI]/executable_ + docExt;
            docFile.expand();

            // Check the existence of the file
            if (isFile(docFile))
            {
                found = true;
                break;
            }
        }
    }

    if (found || httpServer != fileName::null)
    {
        string docBrowser = getEnv("FOAM_DOC_BROWSER");
        if (docBrowser.empty())
        {
            docDict.lookup("docBrowser") >> docBrowser;
        }

        if (found)
        {
            docBrowser += " file://" + docFile;
        }
        else
        {
            docBrowser += " " + httpServer;
        }

        Info<< "Show documentation: " << docBrowser.c_str() << endl;

        system(docBrowser);
    }
    else
    {
        Info<< nl
            << "No documentation found for " << executable_
            << ", but you can use -help to display the usage\n" << endl;
    }
}


bool Foam::argList::check(bool checkArgs, bool checkOpts) const
{
    bool ok = true;

    if (Pstream::master())
    {
        if (checkArgs && args_.size() - 1 != validArgs.size())
        {
            FatalError
                << "Wrong number of arguments, expected " << validArgs.size()
                << " found " << args_.size() - 1 << endl;
            ok = false;
        }

        if (checkOpts)
        {
            forAllConstIter(HashTable<string>, options_, iter)
            {
                if
                (
                    !validOptions.found(iter.key())
                 && !validParOptions.found(iter.key())
                )
                {
                    FatalError
                        << "Invalid option: -" << iter.key() << endl;
                    ok = false;
                }
            }
        }

        if (!ok)
        {
            printUsage();
        }
    }

    return ok;
}


bool Foam::argList::checkRootCase() const
{
    if (!fileHandler().isDir(rootPath()))
    {
        FatalError
            << executable_
            << ": cannot open root directory " << rootPath()
            << endl;

        return false;
    }

    fileName pathDir(fileHandler().filePath(path()));

    if (pathDir.empty() && Pstream::master())
    {
        // Allow slaves on non-existing processor directories, created later
        // (e.g. redistributePar)
        FatalError
            << executable_
            << ": cannot open case directory " << path()
            << endl;

        return false;
    }

    return true;
}


// ************************************************************************* //
