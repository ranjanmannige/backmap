package require Orient
namespace import Orient::orient
# ... load your molecules and make a selection ...
set sel [atomselect top "all"]

# show/calc the principal axes
set I [Orient::calc_principalaxes $sel] 
# rotate axis 0 to match X
set A [orient $sel [lindex $I 2] {1 0 0}]
$sel move $A
# recalc principal axes to check
set I [Orient::calc_principalaxes $sel]  
# rotate axis 1 to match Y
set B [orient $sel [lindex $I 1] {0 1 0}] 
$sel move $B

set nf [molinfo top get numframes]
for {set i 0} {$i < $nf-1} {incr i} { 
	set selfr [atomselect top "all" frame $i]
	$selfr move $A
	$selfr move $B
}
