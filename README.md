# BadgerDb-BTree-Index

In this programming project, we implement a B+Tree index in our database system. A B+Tree is a balanced search tree in which the internal pages direct the search and leaf pages contain the actual data entries. The index provides fast data retrieval without needing to search every row in a database table, enabling rapid random lookups and efficient scans of ordered records. Our implementation support thread-safe search, insertion, deletion (including splitting and merging nodes), and an iterator to support in-order leaf scans.

################################################################################
# BadgerDB quick start guide                                                   #
################################################################################

################################################################################
# Building the source and documentation                                        #
################################################################################

To build the source:
  $ make

To build the real API documentation (requires Doxygen):
  $ make doc

To view the documentation, open docs/index.html in your web browser after
running make doc.

################################################################################
# Prerequisites                                                                #
################################################################################

If you are running this on a CSL instructional machine, these are taken care of.

Otherwise, you need:
 * a modern C++ compiler (gcc version 4.6 or higher, any recent version of clang)
 * doxygen (version 1.4 or higher)
