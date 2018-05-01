%{
/*---------------------------------*- C -*-----------------------------------*\
 =========                   |
 \\      /   F ield          | OpenFOAM: The Open Source CFD Toolbox
  \\    /    O peration      |
   \\  /     A nd            | Copyright (C) 2011-2018 OpenFOAM Foundation
    \\/      M anipulation   |
------------------------------------------------------------------------------
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
    wmkdep

Description
    A fast dependency list generator that emulates the behaviour and the
    output of cpp -M. However, the output contains no duplicates and
    is approx. 40% faster than cpp.

    The algorithm uses flex to scan for includes and searches the files
    found.  Each file is entered into a hash table so that files are scanned
    only once.  This is why this program is faster than cpp.

Usage
    wmkdep [ -Idir ... -Idir ] filename

\*---------------------------------------------------------------------------*/

#define HASH_TABLE_SIZE 500
#define REPLACEMENTS_SIZE 10
#define INITIAL_MAX_N_FILES 1000

#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <dirent.h>
#include <errno.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/* char* entry in hash table */
struct HashEntry
{
    char* name;
    struct HashEntry* next;
};

/* String search/replace pair */
struct searchReplace
{
    char* search;
    size_t searchLen;
    char* replace;
    size_t replaceLen;
};

int nDirectories = 0;
char** directories;
char* sourceFile = NULL;
char* depFilePath = NULL;
char* depFileName = NULL;

int nReplacements = 0;
struct searchReplace replacements[REPLACEMENTS_SIZE];

/* Set of files already visited */
struct HashEntry* visitedFiles[HASH_TABLE_SIZE];

/* List of dependency files */
int nFiles = 0;
int currentFile = 0;
int maxNfiles = INITIAL_MAX_N_FILES;
char** files;

/* Current path of the dependency file */
const char* currentPath = NULL;

void nextFile(const char* fileName);
int lookUp(struct HashEntry** hashTable, const char* p);
void addFile(char* pathName);
void openFile(const char* pathName);
char* strRep(char* str, struct searchReplace* sr);
char* substitutePath(char* pathName);
void printFile(FILE* file, const char* pathName);

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#undef yywrap        /* Sometimes a macro by default */

%}

%x CMNT CFNAME
%%

"//".*\n                  ;               /* Remove c++ style line comments */

"/*"                      BEGIN(CMNT);    /* Start removing c style comment */
<CMNT>.|\n                ;
<CMNT>"*/"                BEGIN(INITIAL); /* End removing c style comment */

^[ \t]*#[ \t]*include[ \t]+\"   BEGIN(CFNAME);  /* c-file name */
<CFNAME>[^"\n ]*        { BEGIN(INITIAL); nextFile(yytext); } /*"*/

.|\t|\n                   ;

%%


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

int main(int argc, char* argv[])
{
    if (argc == 1)
    {
        fprintf(stderr, "input file not supplied\n");
        exit(1);
    }

    sourceFile = strdup(argv[argc-2]);
    depFilePath = strdup(argv[argc-1]);

    char* basePos = strrchr(sourceFile, '/');
    if (basePos == NULL)
    {
        basePos = sourceFile;
    }
    else
    {
        basePos++;
    }

    char* dotPos = strrchr(sourceFile, '.');
    if(dotPos == NULL || dotPos < basePos)
    {
        fprintf
        (
            stderr,
            "cannot find extension in source file name %s\n",
            sourceFile
        );
        exit(1);
    }

    /* Build list of string replacements */
    int i;
    for (i = 1; i < argc - 1; i++)
    {
        if (strncmp(argv[i], "-R", 2) == 0)
        {
            replacements[nReplacements].search = strdup(argv[++i]);
            replacements[nReplacements].searchLen =
                strlen(replacements[nReplacements].search);

            replacements[nReplacements].replace = strdup(argv[++i]);
            replacements[nReplacements].replaceLen =
                strlen(replacements[nReplacements].replace);

            nReplacements++;
        }
    }

    /* Count number of -I directories */
    nDirectories = 0;
    for (i = 1; i < argc - 1; i++)
    {
        if (strncmp(argv[i], "-I", 2) == 0)
        {
            if (strlen(argv[i]) > 2)
            {
                nDirectories++;
            }
        }
    }

    directories = (char**)malloc(sizeof(char*)*nDirectories);

    /* Build list of -I directories */
    nDirectories = 0;
    for (i = 1; i < argc - 1; i++)
    {
        if (strncmp(argv[i], "-I", 2) == 0)
        {
            if (strlen(argv[i]) > 2)
            {
                directories[nDirectories++] = strdup(argv[i] + 2);
            }
        }
    }

    depFileName = substitutePath(strdup(depFilePath));

    /* Initialise storage for the list of the dependencies */
    files = (char**)malloc(sizeof(char*)*maxNfiles);

    openFile(sourceFile);
    yylex();


    /* Write dependencies file */

    FILE* depFile;
    if (!(depFile = fopen(depFilePath, "w")))
    {
        fprintf
        (
            stderr,
            "could not open dependencies file %s "
            "for source file %s due to %s\n",
            depFilePath, sourceFile, strerror(errno)
        );
        exit(1);
    }

    fprintf(depFile, "%s: \\\n", depFileName);
    printFile(depFile, sourceFile);
    for (i = 0; i < nFiles; i++)
    {
        files[i] = substitutePath(files[i]);
        printFile(depFile, files[i]);
    }

    fputs("\n", depFile);

    /* Write the dummy rules for the dependencies */
    for (i = 0; i < nFiles; i++)
    {
        fprintf(depFile, "%s :\n", files[i]);
    }

    fputs("\n", depFile);

    /* Clean-up storage */

    for (i = 0; i < nDirectories; i++)
    {
        free(directories[i]);
    }
    free(directories);

    free(sourceFile);
    free(depFilePath);
    free(depFileName);

    for (i = 0; i < nFiles; i++)
    {
        free(files[i]);
    }
    free(files);

    return 0;
}


/* Add a directory name to the file name */
char* addDirectoryName(const char* dirName, const char* fileName)
{
    char* pathName = (char*)malloc(strlen(dirName) + strlen(fileName) + 2);
    strcpy(pathName, dirName);

    if (dirName[strlen(dirName)-1] != '/')
    {
        strcat(pathName, "/");
    }

    strcat(pathName, fileName);

    return pathName;
}


/* Find path to next file and add it to the list */
void nextFile(const char* fileName)
{
    if (lookUp(visitedFiles, fileName))
    {
        return;   /* Already existed (did not insert) */
    }

    /* Check if the file has same path as the last file read */
    if (currentPath != NULL)
    {
        char* pathName = addDirectoryName(currentPath, fileName);

        if (access(pathName, R_OK ) != -1)
        {
            addFile(pathName);
            return;
        }

        free(pathName);
    }

    if (access(fileName, R_OK ) != -1)
    {
        addFile(strdup(fileName));
        currentPath = NULL;
    }
    else
    {
        int i;
        for (i = 0; i < nDirectories; i++)
        {
            char* pathName = addDirectoryName(directories[i], fileName);

            if (access(pathName, R_OK ) != -1)
            {
                addFile(pathName);
                currentPath = directories[i];

                return;
            }

            free(pathName);
        }

        if (nDirectories == 0)
        {
            fprintf
            (
                stderr,
                "could not open file %s for source file %s\n",
                fileName, sourceFile
            );
        }
        else
        {
            fprintf
            (
                stderr,
                "could not open file %s for source file %s due to %s\n",
                fileName, sourceFile, strerror(errno)
            );
        }

        /* Only report error on the first occurrence */
        lookUp(visitedFiles, fileName);
    }
}


/* Lookup name in hash table.
   If found - return 1
   If not found - insert in table and return 0
*/
int lookUp(struct HashEntry** hashTable, const char* p)
{
    int ii = 0;
    struct HashEntry* n;
    struct HashEntry* nn;

    /* Hash */
    const char* pp = p;
    while (*pp) ii = ii<<1 ^ *pp++;
    if (ii < 0) ii = -ii;
    ii %= HASH_TABLE_SIZE;

    /* Search */
    for (n = hashTable[ii]; n; n = n->next)
    {
        if (strcmp(p, n->name) == 0)
        {
            /* Entry found so return true */
            return 1;
        }
    }

    /* Insert */
    nn = (struct HashEntry*)malloc(sizeof(struct HashEntry));
    nn->name = strdup(p);
    nn->next = hashTable[ii];
    hashTable[ii] = nn;

    /* Entry not found, and therefore added. return false */
    return 0;
}


/* Add file to list */
void addFile(char* pathName)
{
    if (nFiles == maxNfiles - 1)
    {
        maxNfiles *= 2;
        files = (char**)realloc(files, sizeof(char*)*maxNfiles);
    }

    files[nFiles++] = pathName;
}


/* Open file and set yyin */
void openFile(const char* pathName)
{
    if (!(yyin = fopen(pathName, "r")))
    {
        fprintf
        (
            stderr,
            "could not open file %s for source file %s due to %s\n",
            pathName, sourceFile, strerror(errno)
        );
    }
}


/* String search/replace */
char* strRep(char* str, struct searchReplace* sr)
{
    char* searchStart = strstr(str, sr->search);

    if (searchStart  != NULL)
    {
        if (sr->replaceLen > sr->searchLen)
        {
            const size_t start = str - searchStart;
            str = realloc
            (
                str,
                strlen(str) + sr->replaceLen - sr->searchLen + 1
            );
            searchStart = str + start;
        }

        const size_t tailLen = strlen(searchStart + sr->searchLen);

        memmove
        (
            searchStart + sr->replaceLen,
            searchStart + sr->searchLen,
            tailLen + 1
        );

        memcpy(searchStart, sr->replace, sr->replaceLen);
    }

    return str;
}


/* Substitute path components with command-line replacements */
char* substitutePath(char* pathName)
{
    if (nReplacements)
    {
        int i;
        for (i = 0; i < nReplacements; i++)
        {
            pathName = strRep(pathName, &replacements[i]);
        }
    }

    return pathName;
}


/* Print file path to the dependencies file */
void printFile(FILE* file, const char* pathName)
{
    fprintf(file, "%s \\\n", pathName);
}


/* The lexer calls yywrap to handle EOF conditions */
int yywrap()
{
    /* Close the current file which has just reached EOF */
    fclose(yyin);

    if (currentFile == nFiles) /* There are no more files in the list */
    {
        /* Return 1 to inform lex finish now that all files have been read */
        return 1;
    }
    else
    {
        /* Open the next file on the list */
        openFile(files[currentFile++]);

        /* Return to the normal state for the next file */
        BEGIN(INITIAL);

        /* Return 0 to inform lex to continue reading */
        return 0;
    }
}


/*****************************************************************************/
