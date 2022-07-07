#!/usr/bin/env python

import gzip, shutil, sys

print ('\ngzip.py replacement script')
print ('(C) 2022 Christian Wuensch\n')

decompress = False

# Process command line arguments
opts = [opt for opt in sys.argv[1:] if opt.startswith("-")]
args = [arg for arg in sys.argv[1:] if not arg.startswith("-")]

if "-d" in opts:
  decompress = True

if len(args) >= 1:
  infile = args[0]
else:
  print ('Usage: gzip [-cfd] [FILE]...\n')
  print ('Compress FILEs (or stdin)\n')
  print ('\t-d\tDecompress\n')
  exit(1)

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
