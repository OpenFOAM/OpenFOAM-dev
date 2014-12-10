# -----------------------------------------------------------------------------
# Script
#     doxyFilter.sed
#
# Description
#     Transform human-readable tags such as 'Description' into the Doxygen
#     equivalent
# -----------------------------------------------------------------------------

# new FSF address
/^License/,/\*\//{
/^License/,\%http://www.gnu.org/licenses%{
s?^License.*?\*\/\
\/\*! \\file %filePath%\
<b>Original source file</b> <a href="%filePath%">%fileName%</a>\
\
\
\
\
\
\
\
\
?
/^    /d
}


# remove entry
/^Application *$/{
N
N
d
}


# remove entry
/^Global *$/{
N
N
d
}


# Primitive
#     typename
# =>
# \\relates typename
#
/^Primitive *$/,/^[^ ]/{
s/^Primitive *$//
s/^    /\\relates /
}


# Class
#     Foam::className
# =>
# \\class Foam::className
#
/^Class *$/,/^[^ ]/{
s/^Class *$//
s/^    /\\class /
}


# Group
#     groupName
# =>
# \ingroup namespaceName
#
/^Group *$/,/^[^ ]/{
s/^Group//
s/^    /\\ingroup /
}


# Namespace
#     namespaceName
# =>
# \namespace namespaceName
#
/^Namespace *$/,/^[^ ]/{
s/^Namespace//
s/^    /\\namespace /
}


# Typedef
#     Foam::def
# =>
# \typedef Foam::def
/^Typedef *$/,/^[^ ]/{
s/^Typedef//
s/^    /\\typedef /
}


# add anchor and use \brief
# the first paragraph will be 'brief' and the others 'detail'
/^Description *$/,/^[^ ]/{
/^Description/c\
<a class="anchor" name="Description"></a> \\brief
s/^    //
}

/^Usage *$/,/^[^ ]/{
/^Usage/c\
\\par Usage
s/^    //
}


/^See *Also *$/,/^[^ ]/{
/^See *Also/c\
\\see
s/^    //
}

/^Note *$/,/^[^ ]/{
/^Note/c\
\\note
s/^    //
}


# remove ToDo paragraph to avoid them showing on related pages
/^To[Dd]o *$/,/^[^ ]/{
s/^To[Dd]o *$//
s/^    .*//
}


/^Warning *$/,/^[^ ]/{
/^Warning/c\
\\warning
s/^    //
}


/^Deprecated *$/,/^[^ ]/{
/^Deprecated/c\
\\deprecated
s/^    //
}


/^SourceFiles *$/,/^$/{
s?SourceFiles?\\par Source files\
<ul><li><a href="%filePath%">%fileName%</a></li>?
s? *\([a-zA-Z0-9]*\.[a-zA-Z]*\)?  <li><a href="%dirName%/\1">\1</a></li>?
s?^$?</ul>?
}

/fileName%<\/a><\/li>$/{
N
s?\n$?</ul>?g
s/<\/li>\n/<\/li> /
s? *\([a-zA-Z0-9]*\.[a-zA-Z]*\)?  <li><a href="%dirName%/\1">\1</a></li>?
}

s/.*\*\//\*\//


# convert /heading in source files to bold font and add some space
s#\\heading \(.*\)#<br><b>\1</b>#g

# add a linebreak
s#\\linebreak#<br>#g

}

# -----------------------------------------------------------------------------
