prepmate
========

processing mate-pair data from Illumina's nextera protocol.

Inspired by nextclip.

to build simply type

	make

on any unix system

requires that zlib be installed on your system.  Optionally
you can also read from bzip2 formatted files by typing

	make BZIP2=1

when building

prepmate will output five files in gzip format.  Two for
each of the mate-pair and paired-end reads produced 
during processing and a singletons files.

There is also a log file describing the alignments of the 
adaptor sequence to the reads and there eventual fate in 
the dataset.