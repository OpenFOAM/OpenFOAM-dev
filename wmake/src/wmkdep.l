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

#define FILE_STACK_SIZE 300
#define HASH_TABLE_SIZE 500

#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <dirent.h>
#include <errno.h>

void nextFile(const char* fileName);
void importFile(const char* fileName);
void importDir(const char* dirName);

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


/* char* entry in hash table */
struct HashEntry
{
    char* name;
    struct HashEntry* next;
};


/*
 * Lookup name in hash table.
 * If found - return 1
 * If not found - insert in table and return 0
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
    for (n=hashTable[ii]; n; n=n->next)
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


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

int nDirectories = 0;
char** directories;
char* sourceFile = NULL;
char* depFile = NULL;

/* Set of files already visited */
struct HashEntry* visitedFiles[HASH_TABLE_SIZE];

/* Buffer pointer stack counter */
int currentBuffer = 0;

/* Buffer pointer stack */
YY_BUFFER_STATE buffers[FILE_STACK_SIZE];

/* Directory paths for the loaded files */
const char* bufferPaths[FILE_STACK_SIZE];


int main(int argc, char* argv[])
{
    char *basePos, *dotPos;
    int i;

    if (argc == 1)
    {
        fprintf(stderr, "input file not supplied\n");
        exit(1);
    }

    sourceFile = strdup(argv[argc-1]);

    if ((basePos = strrchr(sourceFile, '/')) == NULL)
    {
        basePos = sourceFile;
    }
    else
    {
        basePos++;
    }

    if
    (
        (dotPos = strrchr(sourceFile, '.')) == NULL
     || (dotPos < basePos)
    )
    {
        fprintf
        (
            stderr,
            "cannot find extension in source file name %s\n",
            sourceFile
        );
        exit(1);
    }

    /* Count number of -I directories */
    nDirectories = 0;
    for (i = 1; i < argc; i++)
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
    for (i = 1; i < argc; i++)
    {
        if (strncmp(argv[i], "-I", 2) == 0)
        {
            if (strlen(argv[i]) > 2)
            {
                directories[nDirectories++] = strdup(argv[i] + 2);
            }
        }
    }

    /*
     * Initialise depFile to zero and use strncat rather than strncpy
     * because there is a bug in the SGI strncat that if 0 precedes the '.'
     * it inserts a space
     */
    depFile = (char*)malloc(strlen(sourceFile) + 20);
    depFile[0] = 0;
    strcat(depFile, "$(OBJECTS_DIR)/");
    strcat(depFile, sourceFile);
    strcat(depFile, ".dep");

    char *objectFile = strdup(basePos);
    objectFile[(dotPos - basePos)/sizeof(char)] = 0;

    /* printf("$(OBJECTS_DIR)/%s.o: %s\n", objectFile, depFile); */
    printf("%s: \\\n", depFile);
    free(objectFile);

    /* Initialize buffer path for currentBuffer */
    bufferPaths[currentBuffer] = NULL;

    nextFile(sourceFile);
    yylex();

    puts("\n");

    for (i = 0; i < nDirectories; i++)
    {
        free(directories[i]);
    }
    free(directories);

    free(sourceFile);
    free(depFile);

    return 0;
}


/*
 * Add a directory name to the file name
 */
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


/*
 * Open a file and create buffer and put on stack
 */
void nextFile(const char* fileName)
{
    if (lookUp(visitedFiles, fileName))
    {
        return;   /* Already existed (did not insert) */
    }

    if (currentBuffer >= FILE_STACK_SIZE)
    {
        fprintf
        (
            stderr,
            "depth of file search exceeds stack size %d "
            "while opening %s for file %s\n",
            FILE_STACK_SIZE, fileName, sourceFile
        );
        exit(1);
    }

    /* Pointer to new file which is set if the file is successfully opened */
    FILE* newyyin = NULL;

    /* Check if the file has same path as the file in the current buffer */
    if (bufferPaths[currentBuffer] != NULL)
    {
        char* pathName = addDirectoryName(bufferPaths[currentBuffer], fileName);

        if ((newyyin = fopen(pathName, "r")))
        {
            printf("%s \\\n", pathName);

            buffers[currentBuffer++] = YY_CURRENT_BUFFER;
            bufferPaths[currentBuffer] = bufferPaths[currentBuffer-1];

            yy_switch_to_buffer(yy_create_buffer(newyyin, YY_BUF_SIZE));

            free(pathName);

            return;
        }

        free(pathName);
    }

    if (!(newyyin = fopen(fileName, "r")))
    {
        int d;
        for (d=0; d<nDirectories; d++)
        {
            char* pathName = addDirectoryName(directories[d], fileName);

            if ((newyyin = fopen(pathName, "r")))
            {
                printf("%s \\\n", pathName);

                buffers[currentBuffer++] = YY_CURRENT_BUFFER;
                bufferPaths[currentBuffer] = directories[d];

                yy_switch_to_buffer(yy_create_buffer(newyyin, YY_BUF_SIZE));

                free(pathName);

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

        fflush(stdout);
        fflush(stderr);

        /* Only report the first occurrence */
        lookUp(visitedFiles, fileName);
    }
    else
    {
        printf("%s \\\n", fileName);
        fflush(stdout);

        buffers[currentBuffer++] = YY_CURRENT_BUFFER;
        bufferPaths[currentBuffer] = NULL;

        yy_switch_to_buffer(yy_create_buffer(newyyin, YY_BUF_SIZE));
    }
}


/*
 * Replace all '.' with '/'
 */
void dotToSlash(char* fileName)
{
    int i, len;
    len = strlen(fileName);

    for (i=0; i<len; i++)
    {
        if (fileName[i] == '.') fileName[i] = '/';
    }
}


/*
 * The lexer calls yywrap to handle EOF conditions
 */
int yywrap()
{
    /* Close the file for the buffer which has just reached EOF */
    /* This causes strange problems on some systems
    fclose(yyin);
    yyin = 0;
    */

    /* Delete the buffer */
    yy_delete_buffer(YY_CURRENT_BUFFER);

    /* Set buffer counter to previous buffer */
    currentBuffer--;

    if (currentBuffer >= 0) /* If buffer counter refers to a valid file */
    {
        /* Reset input buffer to the previous buffer on the stack */
        yy_switch_to_buffer(buffers[currentBuffer]);

        /* Return to the normal state for the previous buffer on the stack */
        BEGIN(INITIAL);

        /* Return 0 to inform lex to continue reading */
        return 0;
    }
    else                 /* Else there are no more buffers on the stack */
    {
        /* Return 1 to inform lex finish now that all buffers have been read */
        return 1;
    }
}


/*****************************************************************************/
