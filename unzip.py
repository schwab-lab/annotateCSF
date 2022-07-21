#!/usr/bin/env python

import gzip, shutil, sys

print ('\ngzip.py replacement script')
print ('(C) 2022 Christian Wuensch\n')

# Process command line arguments
opts = [opt for opt in sys.argv[1:] if opt.startswith("-")]
args = [arg for arg in sys.argv[1:] if not arg.startswith("-")]

if len(args) >= 1:
  infile = args[0]
else:
  print ('Usage: python ' + sys.argv[0] + ' [-cd] [FILE]\n')
  print ('Compress FILE\n')
  print ('\t-c\tCompress')
  print ('\t-d\tDecompress\n')
  exit(1)

decompress = (infile.endswith('.gz'))
if "-c" in opts:
  decompress = False
if "-d" in opts:
  decompress = True

if (decompress):
  outfile = infile[0:-3]
  print ('Decompressing "' + infile + '" to "' + outfile +'"...')

  with gzip.open(infile, 'rb') as f_in:
    with open(outfile, 'wb') as f_out:
      shutil.copyfileobj(f_in, f_out)
else:
  outfile = infile + '.gz'
  print ('Compressing "' + infile + '" to "' + outfile +'"...')

  with open(infile, 'rb') as f_in:
    with gzip.open(outfile, 'wb') as f_out:
      shutil.copyfileobj(f_in, f_out)
