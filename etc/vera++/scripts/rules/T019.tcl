#!/usr/bin/tclsh
# control structures should have complete curly-braced block of code

foreach fileName [getSourceFileNames] {

    set state "start"
    set prev ""
    set pp_pragma_line -1
    foreach token [getTokens $fileName 1 0 -1 -1 {for if while do else leftparen rightparen leftbrace rightbrace semicolon pp_pragma}] {
        set type [lindex $token 3]
        set line [lindex $token 1]

        if {$state == "control"} {
            if {$type == "leftparen"} {
                incr parenCount				
            } elseif {$type == "rightparen"} {
                incr parenCount -1
                if {$parenCount == 0} {
                    set state "expectedblock"
                }
            }

        } elseif {$state == "expectedblock"} {
            if {$prev == "else" && $type == "if" } {
              # skip
            } elseif {$type != "leftbrace"} {
                set line [lindex $token 1]
                report $fileName $line "full block {} expected in the control structure"
            }
            set state "block"
        }

        if {$type == "pp_pragma"} {
          set pp_pragma_line $line
        } elseif {$pp_pragma_line != $line} {
            if {$type == "for" || $type == "if"} {
                set parenCount 0
                set state "control"
            } elseif {$type == "do" || $type == "else"} {
                set state "expectedblock"
            } elseif {$type == "while" && $prev != "rightbrace"} {
                set parenCount 0
                set state "control"
            }
        }
        set prev $type
    }
}
