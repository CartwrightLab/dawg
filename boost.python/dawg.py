import dawg_python

akita = dawg_python.DawgWalker(b"../examples/basic-dna.dawg", b"fasta:-", 10, 121)
akita.bark(b"bark bark bark")
akita.run()
