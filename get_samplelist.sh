##########
# Alisa Vershinina
# 12 April 2017
# Goal: get a list of names by a pattern. 
# Usage: 
#	$1 - pattern
#	$2 - d (dirs) or -f (files)
##########

for x in $(find . -name "$1*" -type $2)
do
	echo "$(cut -d _ -f 1,2 <<< $(basename $x))"
done | sort -u | while read y; do echo -n "$y "; done
echo ""
