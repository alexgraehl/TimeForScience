#!/bin/sh

# $1: search for this...
# $2: replace with this...

# $3: search in this file...

#perl -i -p -e 's/$1/$2/g;' $3
sed s/$1/$2/g $3 > $3.tempfile
mv -f $3.tempfile $3

# $4: and replace with this file...

# sed s/$1/$2/g $3 > $4

