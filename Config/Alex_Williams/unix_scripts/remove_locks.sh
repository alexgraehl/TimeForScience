#!/bin/sh

# Delete lock debris on cvs

echo "Removing locks (1)..."
find /projects/compbiousr/cvsroot/svdo/pymol-cvs -name "#cvs.lock*" -print | xargs rm -fdr
echo "Removing locks (2)..."
find /projects/compbiousr/cvsroot/svdo/pymol-cvs -name "#cvs.lock.*" -print | xargs rm -f
echo "Removing locks (3)..."
find /projects/compbiousr/cvsroot/svdo/pymol-cvs -name "#cvs.rfl.*" -print | xargs rm -f
echo "Removing locks (4)..."
find /projects/compbiousr/cvsroot/svdo/pymol-cvs -name "#cvs.wfl.*" -print | xargs rm -f
echo "Removing locks (5)..."
find /projects/compbiousr/cvsroot/svdo/pymol-cvs -name "#cvs.pfl.*" -print | xargs rm -f
echo "Removing locks (6)..."
find /projects/compbiousr/cvsroot/svdo/pymol-cvs -name "#cvs.tfl.*" -print | xargs rm -f
echo "Removing locks (7)..."
find /projects/compbiousr/cvsroot/svdo/pymol-cvs -name "cvsloc" -print | xargs rm -fdr
echo "Removing locks (8)..."

