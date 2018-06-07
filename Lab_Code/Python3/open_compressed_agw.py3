
# In order to use this function, just do:
#  import open_compressed_agw <-- this does not work
# Then you can type 'open_compressed_agw' to transparently
# handle compressed files

import gzip
import bz2
import re


def open_compressed_agw(filename, mode):
	if re.search(r"[.](gz|gzip|Z)$", filename, flags=re.IGNORECASE):
		return gzip.open(filename=filename, mode=mode)
	if re.search(r"[.](bz2|bzip2)$", filename, flags=re.IGNORECASE):
		return bz2.BZ2File(filename=filename, mode=mode)
	else:
		return open(filename=filename, mode=mode)


