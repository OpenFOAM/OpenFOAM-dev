#!/usr/bin/tclsh
# Keywords catch, for, if and while should be followed by a single space or
# a new line

foreach f [getSourceFileNames] {
    set pp_pragma_line -1
    foreach t [getTokens $f 1 0 -1 -1 {catch for if switch while pp_pragma}] {
        set keyword [lindex $t 0]
        set line [lindex $t 1]
        set column [lindex $t 2]
        set type [lindex $t 3]
        if {$type == "pp_pragma"} {
          set pp_pragma_line $line
        } elseif {$pp_pragma_line != $line} {
            set followingTokens [getTokens $f $line [expr $column + [string length $keyword]] [expr $line + 1] -1 {}]
            if {[llength $followingTokens] < 2 & [lindex $followingTokens -1] != "newline"} {
                if {[lindex [lindex $followingTokens 0] 3] != "newline"} {
                    report $f $line "keyword '${keyword}' not followed by a single space or new line"
                }
            } else {
                if {[list [lindex [lindex $followingTokens 0] 0] [lindex [lindex $followingTokens 1] 0]] != [list " " "("]} {
                    report $f $line "keyword '${keyword}' not followed by a single space or new line"
                }
            }
        }
    }
}
